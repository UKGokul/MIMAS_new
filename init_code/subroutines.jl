# src/subroutines.jl
# This file contains utility functions and the main physics subroutines of MIMAS.
# They are called during initialisation and the main time loop. The names and logic
# follow the Fortran original as closely as possible.

# ----------------------------------------------------------------------
# Utility: set_monate – month lengths and leap year
# ----------------------------------------------------------------------
"""
    set_monate(year::Int) -> (Vector{Int}, Int)

Fortran equivalent: `set_monate` (called in main program).
Returns a vector of days in each month (1..12) and the total number of days
in the given year (365 or 366). Leap year rule: year divisible by 4.
(Note: Fortran uses only mod(year,4) == 0, no exception for centuries.)
This is used throughout the code to handle dates and day‑of‑year calculations.
"""
function set_monate(year::Int)
    ntage_monat = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if year % 4 == 0
        ntage_monat[2] = 29
    end
    max_doy = ntage_monat[2] == 29 ? 366 : 365
    return ntage_monat, max_doy
end


# ----------------------------------------------------------------------
# Subroutine: sub_maxmin! – Courant number check and vertical wind capping
# ----------------------------------------------------------------------
"""
    sub_maxmin!(um, vm, wm, S::MIMASState, cfg::Config)

Fortran equivalent: `sub_maxmin` (called after each dynamics read).
Computes Courant numbers for zonal, meridional, and vertical winds:
    ucour_max = max( |um| / dx ) * dttrans
    vcour_max = max( |vm| / dy ) * dttrans
    wcour_max = max( |wm| / dz ) * dttrans
If the vertical Courant number exceeds 0.9, the vertical wind is capped in‑place:
    wm = abs(wm) * 0.9 * dz / (dttrans * wm)
The Courant numbers are printed; if any exceeds 1.0, an alarm message is printed.
"""
function sub_maxmin!(um, vm, wm, S, cfg)
    # Zonal Courant (um on nbneu, igitneu, kgitneu grid)
    umax = 0.0f0
    for k in 1:cfg.kgitneu
        for i in 1:cfg.igitneu
            for j in 1:cfg.nbneu
                val = abs(um[j,i,k] / S.consts.dx[j])
                if val > umax
                    umax = val
                end
            end
        end
    end
    ucour_max = cfg.dttrans * umax

    # Meridional Courant (vm on nbneu1, igitneu, kgitneu grid)
    vmax = 0.0f0
    for k in 1:cfg.kgitneu
        for i in 1:cfg.igitneu
            for j in 1:cfg.nbneu1
                val = abs(vm[j,i,k] / dy)   # dy is constant from constants.jl
                if val > vmax
                    vmax = val
                end
            end
        end
    end
    vcour_max = cfg.dttrans * vmax

    # Vertical Courant (wm on nbneu, igitneu, kgitneu1 grid) and possible capping
    wmax = 0.0f0
    for k in 1:cfg.kgitneu1
        for i in 1:cfg.igitneu
            for j in 1:cfg.nbneu
                wcour = cfg.dttrans * abs(wm[j,i,k]) / dz
                if wcour > 0.9f0
                    # Cap vertical wind as in Fortran
                    wm[j,i,k] = abs(wm[j,i,k]) * 0.9f0 * dz / (cfg.dttrans * wm[j,i,k])
                end
                val = abs(wm[j,i,k] / dz)
                if val > wmax
                    wmax = val
                end
            end
        end
    end
    wcour_max = cfg.dttrans * wmax

    # Print Courant numbers (as in Fortran)
    println("    ucour_max : ", ucour_max)
    println("    vcour_max : ", vcour_max)
    println("    wcour_max : ", wcour_max)
    if max(ucour_max, vcour_max, wcour_max) > 1.0f0
        println("    alarm courant > 1: ", ucour_max, " ", vcour_max, " ", wcour_max)
        # Fortran would optionally stop; we just warn.
    end
end


# ----------------------------------------------------------------------
# Subroutine: init_staub! – particle initialisation (cold start)
# ----------------------------------------------------------------------
"""
    init_staub!(S::MIMASState, cfg::Config)

Fortran equivalent: `sub_init_staub` (called during initialisation when nstart == 0).
Initialises all 40 million particles for a cold start. The distribution is:
  - 10 height bins (k1=1..10) corresponding to 10 different altitude ranges.
  - 25 radius classes (nklasse=1..25) from 0.9 nm to 3.3 nm in steps of 0.1 nm.
Particle counts per height/radius bin are derived from the `nhisto` array
(5 values representing total counts per radius group, spread across heights).
Longitude (xfeld) is uniformly distributed from 1 to igitneu.
Latitude (yfeld) uses a polar‑enhanced distribution: rand in [sin55°, 1] transformed to
   actual latitude via asind, then shifted to start at 55° (index 1 corresponds to ~55°N).
Height (zfeld) depends on the height bin and on pre‑computed indices j_mesopause and
   j_zpress86km. 80% of particles are placed between 85.5–86.5 km, 20% above that,
   with a linear distribution in altitude.
The routine fills S.eis.xfeld, yfeld, zfeld, rinit, and sets rfeld = rinit.
Also stores the latitude‑distribution arrays (bbreite, bsum, sfactor, invfactor) in S
   for later use in particle birth control.
"""
function init_staub!(S, cfg)
    if cfg.nstart == 1
        println("    Restart: skipping particle initialisation.")
        return
    end

    println("    Initialising particles (cold start)...")

        # --- Test mode: if ntrac is small, use a simple uniform distribution ---
    if cfg.ntrac <= 10000
        println("    Using TEST distribution with $(cfg.ntrac) particles.")
        for n in 1:cfg.ntrac
            # Radius between 0.9 and 3.3 nm
            S.eis.rinit[n] = 0.9f0 + 2.4f0 * rand(Float32)
            # Longitude between 1 and igitneu
            S.eis.xfeld[n] = 1.0f0 + rand(Float32) * (cfg.igitneu - 1.0f0)
            # Latitude between 1 and nbneu
            S.eis.yfeld[n] = 1.0f0 + rand(Float32) * (cfg.nbneu - 1.0f0)
            # Height between 1 and kgitneu
            S.eis.zfeld[n] = 1.0f0 + rand(Float32) * (cfg.kgitneu - 1.0f0)
        end
        S.eis.rfeld .= S.eis.rinit
        println("    rfeld first 20: ", S.eis.rfeld[1:min(20, end)])
        println("    rfeld min = ", minimum(S.eis.rfeld), " max = ", maximum(S.eis.rfeld))
        println("    init_staub! (test) done.")
        return
    end

    # Constants from Fortran
    xhisto = Float32[0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                     2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                     3.0, 3.1, 3.2, 3.3]   # 25 radius bin centers (nm)
    nhisto = [864900, 117000, 15800, 2100, 200]   # particle counts per radius group (scaled to 10M)

    # ------------------------------------------------------------------
    # Compute latitude bins (bbreite and bsum) used for latitude selection
    # ------------------------------------------------------------------
    fak = 100.0f0 / 10.36195f0
    S.bsum[1] = 5.46616f0
    S.bbreite[1] = 19.0f0   # corresponds to 55° (since 19+36 = 55)
    for i in 2:35
        S.bbreite[i] = S.bbreite[i-1] + 1.0f0
        wert = fak * cosd(54.5f0 + Float32(i))
        S.bsum[i] = S.bsum[i-1] + wert
    end

    # Pre‑compute sfactor and invfactor (for polar‑enhanced latitude distribution)
    S.sfactor = sind(55.0f0)          # sin(55°)
    S.invfactor = 1.0f0 - S.sfactor   # 1 - sin(55°)
    println("    sfactor: ", S.sfactor, " invfactor: ", S.invfactor)
    # Check: S.invfactor + S.sfactor should be 1.0

    # ------------------------------------------------------------------
    # Build nhunt array (10 heights × 25 radius classes)
    # nhunt(k1, nklasse) gives the number of particles (×10⁷) for that bin.
    # The Fortran code distributes the five nhisto values linearly across the 25 radius classes.
    # We replicate that algorithm exactly.
    # ------------------------------------------------------------------
    nhunt = zeros(Int, 10, 25)
    for k1 in 1:10                # 10 height bins
        for nhaupt in 1:5          # 5 radius groups (each covering 5 radius classes)
            nwerta = nhisto[nhaupt]
            nwertb = (nhaupt < 5) ? nhisto[nhaupt+1] : 1

            wert1 = Float32(nwertb) + 5.0f0 * (nwerta - nwertb)
            wert2 = Float32(nwertb) + 4.0f0 * (nwerta - nwertb)
            wert3 = Float32(nwertb) + 3.0f0 * (nwerta - nwertb)
            wert4 = Float32(nwertb) + 2.0f0 * (nwerta - nwertb)
            wert5 = Float32(nwertb) + 1.0f0 * (nwerta - nwertb)

            sumw = wert1 + wert2 + wert3 + wert4 + wert5
            scalewert = nwerta / sumw

            nwert1 = round(Int, wert1 * scalewert)
            nwert2 = round(Int, wert2 * scalewert)
            nwert3 = round(Int, wert3 * scalewert)
            nwert4 = round(Int, wert4 * scalewert)
            nwert5 = round(Int, wert5 * scalewert)
            nsum = nwert1 + nwert2 + nwert3 + nwert4 + nwert5
            nwert1 += nwerta - nsum   # adjust to exact total

            idx = (nhaupt-1)*5 + 1
            nhunt[k1, idx]   = nwert1
            nhunt[k1, idx+1] = nwert2
            nhunt[k1, idx+2] = nwert3
            nhunt[k1, idx+3] = nwert4
            nhunt[k1, idx+4] = nwert5
        end
    end

    # Print a small sample of nhunt to verify (first height bin, first few radius classes)
    println("    nhunt[1,1:5] = ", nhunt[1,1:5])

    # ------------------------------------------------------------------
    # Fill particles
    # ------------------------------------------------------------------
    nsum = 0
    for k1 in 1:10                # height bins
        for nklasse in 1:25       # radius classes
            nzahl = nhunt[k1, nklasse] * 4   # scale from 10M to 40M
            for _ in 1:nzahl
                nsum += 1
                if nsum > cfg.ntrac
                    error("Too many particles: nsum > ntrac")
                end

                # Radius: rinit = xhisto[nklasse] + 0.1 * rand()
                rando = rand(Float32)
                S.eis.rinit[nsum] = xhisto[nklasse] + 0.1f0 * rando

                # Longitude: xfeld uniformly distributed 1..igitneu
                rando = rand(Float32)
                S.eis.xfeld[nsum] = 1.0f0 + rando * cfg.igitneu
                # Ensure within [1, igitneu] (tiny adjustment to avoid floating‑point boundaries)
                if S.eis.xfeld[nsum] > cfg.igitneu
                    S.eis.xfeld[nsum] = cfg.igitneu - Base.eps(Float32)
                elseif S.eis.xfeld[nsum] < 1.0f0
                    S.eis.xfeld[nsum] = 1.0f0 + Base.eps(Float32)
                end

                # Latitude: polar‑enhanced distribution
                rando = rand(Float32)
                rando = S.sfactor + S.invfactor * rando   # in [sin55°, 1]
                wert = asind(rando)                        # latitude in degrees
                S.eis.yfeld[nsum] = wert - 36.0f0          # shift so that index 1 corresponds to 55°N
                # Clamp to [1, nbneu]
                if S.eis.yfeld[nsum] > cfg.nbneu
                    S.eis.yfeld[nsum] = cfg.nbneu - Base.eps(Float32)
                elseif S.eis.yfeld[nsum] < 1.0f0
                    S.eis.yfeld[nsum] = 1.0f0 + Base.eps(Float32)
                end

                # Height: depends on height bin k1 and the mesopause/86km indices
                jj = round(Int, S.eis.yfeld[nsum])   # latitude index
                jj = clamp(jj, 1, cfg.nbneu)
                rando = rand(Float32)
                if rando > 0.8f0
                    # Upper tail (20%) – height above 86.5 km
                    rando2 = (rando - 0.8f0) * 5.0f0   # rescale to [0,1]
                    kk1 = S.j_mesopause[jj] - (5 + S.j_zpress86km[jj])
                    kk2 = 15
                    if kk1 >= 1
                        kk2 += kk1
                    end
                    S.eis.zfeld[nsum] = Float32(S.j_zpress86km[jj]) - 3.0f0 + Float32(kk2) * rando2
                else
                    # Main 80%: between 85.5 and 86.5 km
                    rando2 = rando * (1.0f0 / 0.8f0)   # rescale to [0,1]
                    S.eis.zfeld[nsum] = Float32(S.j_zpress86km[jj]) - 5.0f0 + 10.0f0 * rando2
                end

                # Clamp zfeld to [1, kgitneu]
                if S.eis.zfeld[nsum] > cfg.kgitneu
                    S.eis.zfeld[nsum] = cfg.kgitneu - Base.eps(Float32)
                elseif S.eis.zfeld[nsum] < 1.0f0
                    S.eis.zfeld[nsum] = 1.0f0 + Base.eps(Float32)
                end
            end
        end
    end

    # Set rfeld = rinit (all particles start as dust, no ice yet)
    S.eis.rfeld .= S.eis.rinit

    # Print summary statistics
    println("    rfeld first 20: ", S.eis.rfeld[1:20])
    println("    rfeld last 20: ", S.eis.rfeld[end-19:end])
    println("    rfeld min = ", minimum(S.eis.rfeld), " max = ", maximum(S.eis.rfeld))
    println("    nsum = ", nsum, " (should be ", cfg.ntrac, ")")
    println("    init_staub! done.")
end


# ----------------------------------------------------------------------
# Interpolation routines (from Fortran sub_intpol, sub_yspline, sub_raspl1)
# ----------------------------------------------------------------------
"""
    sub_intpol(x, y, xi) -> yi

Perform rational spline interpolation. Given arrays x (n) and y (n), compute
interpolated values yi at points xi (ni). Uses sub_raspl1 and sub_yspline.
"""
function sub_intpol(x, y, xi)
    n = length(x)
    ni = length(xi)
    p = fill(4.1f0, 380)
    q = fill(4.1f0, 380)
    y1 = zeros(Float32, 380)
    a = zeros(Float32, 380)
    b = zeros(Float32, 380)
    c = zeros(Float32, 380)
    d = zeros(Float32, 380)

    y1[1] = (y[2] - y[1]) / (x[2] - x[1])
    y1[n] = (y[n] - y[n-1]) / (x[n] - x[n-1])

    sub_raspl1(x, y, p, q, y1, a, b, c, d, n)

    yi = zeros(Float32, ni)
    for i in 1:ni
        yi[i] = sub_yspline(x, xi[i], a, b, c, d, p, q, n)
    end
    return yi
end

"""
    sub_yspline(x, xwert, a, b, c, d, p, q, n) -> ywert

Evaluate the rational spline at a single point.
"""
function sub_yspline(x, xwert, a, b, c, d, p, q, n)
    i = 1
    while i < n && xwert > x[i+1]
        i += 1
    end
    if i == n
        i = n-1
    end
    t = (xwert - x[i]) / (x[i+1] - x[i])
    u = 1.0f0 - t
    ywert = a[i] * u + b[i] * t + c[i] * u^3 / (p[i] * t + 1.0f0) +
            d[i] * t^3 / (q[i] * u + 1.0f0)
    return ywert
end

"""
    sub_raspl1(x, y, p, q, y1, a, b, c, d, n)

Compute spline coefficients.
"""
function sub_raspl1(x, y, p, q, y1, a, b, c, d, n)
    n1 = n - 1
    c[1] = 0.0f0
    d[1] = 0.0f0
    j1 = 1
    p21 = 0.0f0
    qq1 = 0.0f0
    h1 = 0.0f0
    r1 = 0.0f0
    for k in 1:n1
        j2 = k + 1
        pp = p[k]
        qq = q[k]
        pp2 = pp * (pp + 3.0f0) + 3.0f0
        qq2 = qq * (qq + 3.0f0) + 3.0f0
        p22 = 2.0f0 + pp
        q22 = 2.0f0 + qq
        a[k] = x[j2] - x[k]
        h = 1.0f0 / a[k]
        b[k] = 1.0f0 / (p22 * q22 - 1.0f0)
        h2 = h * b[k]
        r2 = h * h2 * (y[j2] - y[k])
        if k != 1
            hq = h1 * qq1
            hp = h2 * p22
            z = 1.0f0 / (hq * (p21 - c[j1]) + hp * q22)
            c[k] = z * hp
            h_val = r1 * qq1 * (1.0f0 + p21) + r2 * pp2 * (1.0f0 + q22)
            if k == 2
                h_val = h_val - hq * y1[1]
            end
            if k == n1
                h_val = h_val - hp * y1[n]
            end
            d[k] = z * (h_val - hq * d[j1])
        end
        j1 = k
        p21 = p22
        qq1 = qq2
        h1 = h2
        r1 = r2
    end
    y1[n1] = d[n1]
    if n1 > 2
        for j1 in 2:n1-1
            k = n - j1
            y1[k] = d[k] - c[k] * y1[k+1]
        end
    end
    for k in 1:n1
        j2 = k + 1
        h = b[k] * (y[j2] - y[k])
        z = b[k] * a[k]
        p2 = 2.0f0 + p[k]
        q2 = 2.0f0 + q[k]
        c[k] = (1.0f0 + q2) * h - z * (y1[j2] + q2 * y1[k])
        d[k] = -(1.0f0 + p2) * h + z * (p2 * y1[j2] + y1[k])
        a[k] = y[k] - c[k]
        b[k] = y[j2] - d[k]
    end
    nothing
end

# ----------------------------------------------------------------------
# H2O initialisation: interpolate from pressure levels to model grid
# ----------------------------------------------------------------------
"""
    h2oinit_zpr_to_zgeo!(S::MIMASState, cfg::Config)

Fortran equivalent: `sub_h2oinit_zpr_to_zgeo`.
For each (j,i) column, interpolate the H2O mixing ratio from the pressure‑level
data (S.h2o.zpress_init, S.h2o.zpress_init_h2o) onto the model's pressure altitude
grid (S.uvwtd48.zm1). The result is stored in S.h2o.hminit3d.
This routine is called at the beginning of each hour.
"""
function h2oinit_zpr_to_zgeo!(S, cfg)
    nbneu = cfg.nbneu
    igitneu = cfg.igitneu
    kgitneu = cfg.kgitneu
    kgitzpr = cfg.kgitzpr

    zpress_init = S.h2o.zpress_init
    zpress_init_h2o = S.h2o.zpress_init_h2o
    zm = S.uvwtd48.zm1   # current pressure altitude

    for j in 1:nbneu
        for i in 1:igitneu
            # Extract the vertical profile of pressure altitude at this (j,i)
            zpressakt = [zm[j, i, k] for k in 1:kgitneu]
            # Interpolate H2O from pressure levels to these altitudes
            fneu = sub_intpol(zpress_init, zpress_init_h2o, zpressakt)
            # Store in hminit3d
            for k in 1:kgitneu
                S.h2o.hminit3d[j, i, k] = fneu[k]
            end
        end
    end
    println("    h2oinit_zpr_to_zgeo! done.")
end


"""
    sub_zenit!(ntg, itag_im_jahr, i_hemi, S, cfg)

Fortran equivalent: `sub_zenit` (called at the start of each hour with ntg=1, and then each timestep).
Computes the zenith angle factor `secchi(j,i)` for each grid point using the solar declination
(oblinoa) and the current time of day. The result is stored in S.h2o.secchi.
Uses the Chapman function approximations as in the original code.
"""
function sub_zenit!(ntg, itag_im_jahr, i_hemi, S, cfg)
    # Local constants (from Fortran data statements)
    skh = 7.0f0
    zgeo_86 = 86.0f0
    za = 1.0606963f0
    zb = 0.55643831f0
    zd = 1.7245609f0
    zccc = 1.0619896f0
    zf = 0.56498823f0
    zg = 0.06651874f0

    # Compute sunris_86 (sunrise angle at 86 km)
    sunris_86 = acos(6373.0f0 / (6373.0f0 + zgeo_86)) * 180.0f0 / pi + 90.0f0

    # Adjust day of year for hemisphere
    itag = itag_im_jahr
    if i_hemi == 1
        itag = mod(itag_im_jahr + 182, 365)
        if itag == 0
            itag = 365
        end
    end

    # Solar declination (from oblinoa array, which is in degrees)
    dek = S.h2o.oblinoa[itag] * pi / 180.0f0

    # Timestep parameters
    dayl = 86400.0f0 / cfg.dttrans   # timesteps per day (960)
    day16 = dayl / cfg.igitneu       # timesteps per longitude bin
    rncom = Float32(ntg)
    zeit = mod(rncom, dayl)

    # Loop over all grid points
    for i in 1:cfg.igitneu
        for j in 1:cfg.nbneu
            # Latitude in radians
            breite = (37.5f0 + (j-1)) * pi / 180.0f0

            # Local time (in timesteps) for this longitude
            ortzeit = zeit + (i-1) * day16
            ortzeit = mod(ortzeit, dayl)
            stuwi = ortzeit / dayl * 2.0f0 * pi   # hour angle

            # Cosine of zenith angle
            coschi = -cos(breite) * cos(dek) * cos(stuwi) + sin(breite) * sin(dek)
            coschi = clamp(coschi, -1.0f0, 1.0f0)

            chi = acos(coschi) * 180.0f0 / pi   # zenith angle in degrees

            # Compute secchi (air mass factor) using Chapman function approximations
            sec = 0.0f0
            if chi <= 75.0f0
                sec = 1.0f0 / coschi
            elseif chi > 75.0f0 && chi < 90.0f0
                x = (6373.0f0 + zgeo_86) / skh
                y = sqrt(0.5f0 * x) * abs(cos(chi * pi / 180.0f0))
                erfy = (za + zb * y) / (zccc + zd * y + y * y)
                if y > 8.0f0
                    erfy = zf / (zg + y)
                end
                sec = sqrt(pi / 2.0f0 * x) * erfy
            elseif chi > sunris_86
                sec = 0.0f0
            else
                x = (6373.0f0 + zgeo_86) / skh
                y = sqrt(0.5f0 * x) * abs(cos(chi * pi / 180.0f0))
                erfy = (za + zb * y) / (zccc + zd * y + y * y)
                if y > 8.0f0
                    erfy = zf / (zg + y)
                end
                arg = x * (1.0f0 - sin(chi * pi / 180.0f0))
                sec = sqrt(2.0f0 * pi * x) * (sqrt(sin(chi * pi / 180.0f0)) * exp(arg) - 0.5f0 * erfy)
            end
            S.h2o.secchi[j,i] = sec
        end
    end
    nothing
end


# ----------------------------------------------------------------------
# Photolysis and nudging (sub_photolyse)
# ----------------------------------------------------------------------
"""
    sub_photolyse!(mnudge, delta_t, itag1, ihemi, xlyman_obs, S, cfg)

Fortran equivalent: `sub_photolyse`.
Computes photodissociation of H2O by Lyman‑α radiation. Updates `S.h2o.hm` in‑place.
If `mnudge == 1`, also applies a nudging term toward the zonal mean `S.h2o.hminit`.
"""

const _PHOTOLYSE_PRINTED = Ref(false)
function sub_photolyse!(mnudge, delta_t, itag1, ihemi, xlyman_obs, S, cfg)
    # ------------------------------------------------------------------
    # 1. Adjust day index for hemisphere
    # ------------------------------------------------------------------
    itag = itag1
    if ihemi == 1
        itag = mod(itag1 + 182, 365)
        if itag == 0
            itag = 365
        end
    end

    # ------------------------------------------------------------------
    # 2. Top‑of‑atmosphere Lyman‑α photon flux (cm⁻² s⁻¹)
    # ------------------------------------------------------------------
    # phitop_ly = xlyman_obs[itag] * 1e11f0 # faster, but f0 causing error in my julia env.
    phitop_ly = xlyman_obs[itag] * Float32(1e11)
    # ------------------------------------------------------------------
    # 3. Constants for reduction factor (Chabrillat & Kockarts 1997)
    # ------------------------------------------------------------------
    b1_h2o = 0.68431f0
    c1_h2o = 8.2214f-21
    b2_h2o = 0.229841f0
    c2_h2o = 1.77556f-20
    b3_h2o = 0.0865412f0
    c3_h2o = 8.22112f-21

    # ------------------------------------------------------------------
    # 4. Compute O2 column (so2neu) – used in reduction factor
    # ------------------------------------------------------------------
    so2neu = zeros(Float32, cfg.nbneu, cfg.igitneu)
    # Start with the top level (k = kgitneu)
    for i in 1:cfg.igitneu
        for j in 13:cfg.nbneu
            xnluft = Float32(1e-6) * S.uvwtd.dm[j, i, cfg.kgitneu] /
                     (xkbolz * S.uvwtd.tm[j, i, cfg.kgitneu])
            xno2 = 0.2f0 * xnluft
            so2neu[j, i] = xno2 * 4.5f0 * Float32(1e5)   # 4.5 km scale height in cm
        end
    end

    # ------------------------------------------------------------------
    # 5. Main loop from top to bottom (k = kgitneu down to 3)
    # ------------------------------------------------------------------
    for k in cfg.kgitneu:-1:3
        for i in 1:cfg.igitneu
            for j in 1:cfg.nbneu
                # Number density of air (cm⁻³) and H2O (cm⁻³)
                xnluft = Float32(1e-6) * S.uvwtd.dm[j, i, k] / (xkbolz * S.uvwtd.tm[j, i, k])
                xnh2o = S.h2o.hm[j, i, k] * Float32(1e-6) * xnluft

                # Update O2 column for this layer (100 m → 1e4 cm)
                xno2 = 0.2f0 * xnluft
                so2neu[j, i] += xno2 * Float32(1e4)

                # Only if sun is up (secchi > 0.1)
                if S.h2o.secchi[j, i] > 0.1f0
                    wert = so2neu[j, i] * S.h2o.secchi[j, i]   # O2 column along slant path

                    # Exponentials (prevent underflow)
                    exply1 = -c1_h2o * wert
                    exply1 = max(exply1, -70.0f0)
                    exply2 = -c2_h2o * wert
                    exply2 = max(exply2, -70.0f0)
                    exply3 = -c3_h2o * wert
                    exply3 = max(exply3, -70.0f0)

                    reduct_h2o = b1_h2o * exp(exply1) +
                                 b2_h2o * exp(exply2) +
                                 b3_h2o * exp(exply3)

                    photodiss_rate_h2o = phitop_ly * photodiss_cross_h2o * reduct_h2o
                    gesamt_photo_h2o = photodiss_rate_h2o * xnh2o * delta_t

                    xnh2o -= gesamt_photo_h2o
                    # Convert back to ppmv
                    S.h2o.hm[j, i, k] = xnh2o * Float32(1e6) / xnluft
                end
            end
        end
    end

    # ------------------------------------------------------------------
    # 6. Nudging toward zonal mean (only if mnudge == 1)
    # ------------------------------------------------------------------
    if mnudge == 1
        xnudge = 90.0f0 / 86400.0f0   # ~0.00104167

        # bfak (latitude‑dependent factor) – defined for j=1..20, else 0.1
        bfak = zeros(Float32, cfg.nbneu)
        for j in 1:cfg.nbneu
            bfak[j] = 0.1f0
        end
        bfak[1] = 1.5f0;  bfak[2] = 1.5f0;  bfak[3] = 1.5f0;  bfak[4] = 1.5f0;  bfak[5] = 1.5f0
        bfak[6] = 1.4f0;  bfak[7] = 1.3f0;  bfak[8] = 1.2f0;  bfak[9] = 1.1f0;  bfak[10] = 1.0f0
        bfak[11] = 1.0f0; bfak[12] = 1.0f0; bfak[13] = 0.9f0; bfak[14] = 0.8f0; bfak[15] = 0.7f0
        bfak[16] = 0.6f0; bfak[17] = 0.5f0; bfak[18] = 0.4f0; bfak[19] = 0.3f0; bfak[20] = 0.2f0

        # hfak (height‑dependent factor) – computed for all k
        hfak = zeros(Float32, cfg.kgitneu)
        for k in 1:5
            hfak[k] = 40.0f0
        end
        for k in 6:10
            hfak[k] = 40.0f0
        end
        hfak[11] = 38.0f0; hfak[12] = 36.0f0; hfak[13] = 34.0f0; hfak[14] = 32.0f0; hfak[15] = 30.0f0
        hfak[16] = 28.0f0; hfak[17] = 26.0f0; hfak[18] = 24.0f0; hfak[19] = 22.0f0; hfak[20] = 20.0f0
        hfak[21] = 18.0f0; hfak[22] = 16.0f0; hfak[23] = 14.0f0; hfak[24] = 12.0f0; hfak[25] = 10.0f0
        hfak[26] =  8.0f0; hfak[27] =  6.0f0; hfak[28] =  4.0f0; hfak[29] =  2.0f0
        for k in 30:cfg.kgitneu
            hfak[k] = 1.0f0
        end
        # Mirror for upper levels? Fortran does: do k=1,29; hfak(kgitneu+1-k)=hfak(k)*0.5; enddo
        # That is for kgitneu = 163, it sets hfak[163- k +1] = hfak[k]*0.5 for k=1..29.
        # So we need to implement that.
        for k in 1:29
            if cfg.kgitneu + 1 - k <= cfg.kgitneu
                hfak[cfg.kgitneu + 1 - k] = hfak[k] * 0.5f0
            end
        end

        # Apply nudging
        for k in 1:cfg.kgitneu
            for j in 1:cfg.nbneu
                faknudge = min(4.0f0, 5.0f0 * bfak[j] * hfak[k])
                for i in 1:cfg.igitneu
                    S.h2o.hm[j, i, k] += xnudge * faknudge * (S.h2o.hminit[j, k] - S.h2o.hm[j, i, k])
                end
            end
        end
    end

    # Diagnostic print (only once, after first nudging call)
    if mnudge == 1 && !_PHOTOLYSE_PRINTED[]
        println("    photolyse (first nudging call):")
        println("      sample hm[26,1,100] = ", S.h2o.hm[26,1,100])
        println("      min hm = ", minimum(S.h2o.hm), " max = ", maximum(S.h2o.hm))
        _PHOTOLYSE_PRINTED[] = true
    end

    nothing
end


# ----------------------------------------------------------------------
# Particle microphysics and transport (sub_tracertransp)
# ----------------------------------------------------------------------
"""
    sub_tracertransp!(S::MIMASState, cfg::Config)

Fortran equivalent: `sub_tracertransp`.
Computes ice particle growth/sublimation, nucleation, and transport.
Updates:
  - S.eis.rfeld (particle radii)
  - S.eis.xfeld, yfeld, zfeld (particle positions)
  - S.h2o.hm (H2O mixing ratio)
Uses lookup tables from S.tab and dynamics from S.uvwtd.
This is the core microphysics routine, called each timestep.
"""
function sub_tracertransp!(S, cfg)
    # ------------------------------------------------------------------
    # Local constants (converted to Float32)
    # ------------------------------------------------------------------
    eps = 1.0f-4
    aconst1 = 1.0f-27 * rhoice * 4.0f0 * pi / (3.0f0 * xmh2o)

    # Pre‑compute advection factors (dttrans_2 / spacing)
    ax = @. cfg.dttrans_2 / S.consts.dx          # vector of length nbneu
    ay = Float32(cfg.dttrans_2) / dy
    az = Float32(cfg.dttrans_2) / dz

    # work1 = 1e-6 * hm * dm   (used in supersaturation calculation)
    S.work1 .= Float32(1e-6) .* S.h2o.hm .* S.uvwtd.dm

    # Zero accumulation arrays
    fill!(S.change, 0.0f0)
    fill!(S.change_n, 0.0f0)

    # ------------------------------------------------------------------
    # First particle loop: microphysics (growth / nucleation)
    # ------------------------------------------------------------------
    # This loop is over all 40 million particles. For performance,
    # we use @inbounds and keep the loop simple. Later we can add
    # threading with Threads.@threads.
    for n in 1:cfg.ntrac
        # Particle indices (integer positions)
        i = floor(Int, S.eis.xfeld[n])
        j = floor(Int, S.eis.yfeld[n])
        k = floor(Int, S.eis.zfeld[n])

        # Ensure indices are within valid range (should be, but safety)
        i = clamp(i, 1, cfg.igitneu)
        j = clamp(j, 1, cfg.nbneu)
        k = clamp(k, 1, cfg.kgitneu)

        rwert0 = S.eis.rfeld[n]
        rinit_n = S.eis.rinit[n]

        tback = S.uvwtd.tm[j, i, k]
        tp = tback

        if rwert0 > rinit_n
            # Ice particle: growth or sublimation
            work1akt = S.work1[j, i, k]

            # Particle temperature adjustment (worka_tab)
            irad = floor(Int, rwert0 + 0.1f0)
            irad = clamp(irad, 1, size(S.tab.worka_tab, 1))
            tp = S.tab.worka_tab[irad, k] * 1.0f0 + tp
            if tp <= 100.0f0
                tp = 100.1f0
            end

            itemp = floor(Int, tp - 99.0f0)
            itemp = clamp(itemp, 1, 71)

            # Kelvin‑corrected saturation vapour pressure
            irad_kelvin = min(300, floor(Int, 10.0f0 * rwert0))
            irad_kelvin = max(1, irad_kelvin)
            vpsatt_tp = S.tab.vpsatt_tab[itemp] *
                        S.tab.faktkelvin_tab[irad_kelvin, itemp]
            satt_tp = work1akt / vpsatt_tp

            # Radius change rate
            drdt = (satt_tp - 1.0f0) * vpsatt_tp * S.tab.sqrt1t_tab[itemp]
            rwert1 = rwert0 + drdt * Float32(cfg.dttrans_2)
            rwert1 = max(rinit_n, rwert1)
            S.eis.rfeld[n] = rwert1

            # Water mass exchanged (in molecules? then converted to concentration)
            backgr = aconst1 * (rwert1^3 - rwert0^3)
            volufak = S.tab.horiwicht_tab[j]
            S.change_n[n] = backgr * volufak * 0.7f0 * 8.0f0 / 10.0f0

        else
            # Dust particle: possible nucleation
            if tback < 155.0f0
                work1akt = S.work1[j, i, k]

                itemp = floor(Int, tp - 99.0f0)
                itemp = clamp(itemp, 1, 71)

                vpsatt_tp = S.tab.vpsatt_tab[itemp]
                irad_kelvin = floor(Int, 10.0f0 * rwert0)
                irad_kelvin = clamp(irad_kelvin, 1, size(S.tab.faktkelvin_tab, 1))
                vpsatt_tp *= S.tab.faktkelvin_tab[irad_kelvin, itemp]

                satt_tp = work1akt / vpsatt_tp

                if satt_tp > 1.0f0
                    drdt = (satt_tp - 1.0f0) * vpsatt_tp * S.tab.sqrt1t_tab[itemp]
                    rwert1 = rwert0 + drdt * Float32(cfg.dttrans_2)
                    rwert1 = max(rinit_n, rwert1)
                    S.eis.rfeld[n] = rwert1

                    backgr = aconst1 * (rwert1^3 - rwert0^3)
                    volufak = S.tab.horiwicht_tab[j]
                    S.change_n[n] = backgr * volufak * 0.7f0 * 8.0f0 / 10.0f0
                end
            end
        end
    end

    # ------------------------------------------------------------------
    # Random numbers for vertical turbulent displacement
    # ------------------------------------------------------------------
    # Generate two random numbers to define the per‑particle sequence
    # (as in Fortran). We'll use simple rand() for now.
    intel_no = rand(1:50)
    intel_offset = rand(0:999999)

    # ------------------------------------------------------------------
    # Second particle loop: advection (with random walk and sedimentation)
    # ------------------------------------------------------------------
    for n in 1:cfg.ntrac
        i = floor(Int, S.eis.xfeld[n])
        j = floor(Int, S.eis.yfeld[n])
        k = floor(Int, S.eis.zfeld[n])

        i = clamp(i, 1, cfg.igitneu)
        j = clamp(j, 1, cfg.nbneu)
        k = clamp(k, 1, cfg.kgitneu)

        # Temperature index for sedimentation
        itemp = floor(Int, S.uvwtd.tm[j, i, k] - 99.0f0)
        itemp = clamp(itemp, 1, 71)

        # Sedimentation velocity (positive downward)
        wfall = S.eis.rfeld[n] * S.tab.sqrtsedi_tab[itemp] / S.uvwtd.dm[j, i, k]

        # Random vertical velocity from turbulence
        # Using the same scheme as Fortran: 0.5*(2*rand -1) * wrichtig_tab
        wrando = 0.5f0 * (2.0f0 * rand(Float32) - 1.0f0) * S.tab.wrichtig_tab[k]

        # Advect in x (zonal)
        S.eis.xfeld[n] += ax[j] * S.uvwtd.um[j, i, k]

        # Advect in y (meridional) – note: vm is on staggered grid, but we use the value at the particle's j.
        S.eis.yfeld[n] += ay * S.uvwtd.vm[j, i, k]

        # Advect in z (vertical) with random walk and sedimentation
        S.eis.zfeld[n] += az * (S.uvwtd.wm[j, i, k] + wrando - wfall)

        # ------------------------------------------------------------------
        # Boundary handling (wrap in x, reflect in y, clamp in z)
        # ------------------------------------------------------------------
        jj = floor(Int, S.eis.yfeld[n])
        if jj > cfg.nbneu
            # Reflect at north pole
            S.eis.yfeld[n] = 2.0f0 * (cfg.nbneu + 1) - S.eis.yfeld[n]
            jj = floor(Int, S.eis.yfeld[n])
            if jj > cfg.nbneu
                S.eis.yfeld[n] -= eps
            end
            # Wrap longitude (x) by half the globe (?) Fortran: xfeld = 1. + amod(xfeld-1.+igitneu/2., igitneu)
            S.eis.xfeld[n] = 1.0f0 + mod(S.eis.xfeld[n] - 1.0f0 + cfg.igitneu/2.0f0, cfg.igitneu)
        elseif jj < 1
            S.eis.yfeld[n] = 1.0001f0
        end

        ii = floor(Int, S.eis.xfeld[n])
        if ii > cfg.igitneu
            S.eis.xfeld[n] -= cfg.igitneu
            ii = floor(Int, S.eis.xfeld[n])
            if ii < 1
                S.eis.xfeld[n] += eps
            end
        elseif ii < 1
            S.eis.xfeld[n] += cfg.igitneu
            ii = floor(Int, S.eis.xfeld[n])
            if ii > cfg.igitneu
                S.eis.xfeld[n] -= eps
            end
        end

        kk = floor(Int, S.eis.zfeld[n])
        if kk >= cfg.kgitneu
            S.eis.zfeld[n] = cfg.kgitneu - eps   # just below top boundary
        elseif kk <= 1
            S.eis.zfeld[n] = 1.0f0 + eps        # just above bottom
        end
    end

    # ------------------------------------------------------------------
    # Accumulate H2O change from particles into the grid
    # ------------------------------------------------------------------
    for n in 1:cfg.ntrac
        i = floor(Int, S.eis.xfeld[n])
        j = floor(Int, S.eis.yfeld[n])
        k = floor(Int, S.eis.zfeld[n])

        i = clamp(i, 1, cfg.igitneu)
        j = clamp(j, 1, cfg.nbneu)
        k = clamp(k, 1, cfg.kgitneu)

        S.change[j, i, k] += S.change_n[n]
    end

    # ------------------------------------------------------------------
    # Zonal mean smoothing at the pole (j = nbneu)
    # ------------------------------------------------------------------
    for k in 17:cfg.kgitneu-3
        dsum = 0.0f0
        for i in 1:cfg.igitneu
            dsum += S.change[cfg.nbneu, i, k]
        end
        dsum /= cfg.igitneu
        for i in 1:cfg.igitneu
            S.change[cfg.nbneu, i, k] = dsum
        end
    end

    # ------------------------------------------------------------------
    # Update H2O mixing ratio
    # ------------------------------------------------------------------
    # hm_new = max(0.01, hm - change * xkbolz * tm / dm)
    @. S.h2o.hm = max(0.01f0, S.h2o.hm - S.change * xkbolz * S.uvwtd.tm / S.uvwtd.dm)

    # ------------------------------------------------------------------
    # Zonal mean smoothing at the pole for hm
    # ------------------------------------------------------------------
    for k in 17:cfg.kgitneu-3
        dsum = 0.0f0
        for i in 1:cfg.igitneu
            dsum += S.h2o.hm[cfg.nbneu, i, k]
        end
        dsum /= cfg.igitneu
        for i in 1:cfg.igitneu
            S.h2o.hm[cfg.nbneu, i, k] = dsum
        end
    end

    # ------------------------------------------------------------------
    # Horizontal smoothing of hm near the pole (9‑point weighted average)
    # ------------------------------------------------------------------
    # Save current hm into work1
    S.work1 .= S.h2o.hm

    for k in 17:cfg.kgitneu-3
        for i in 1:cfg.igitneu
            # Neighbour indices in longitude (periodic)
            i0 = i == 1 ? cfg.igitneu : i-1
            i2 = i == cfg.igitneu ? 1 : i+1

            # j = nbneu (pole)
            j = cfg.nbneu
            j0 = j-1
            j2 = j   # j2 = nbneu (since j+1 would be out of bounds)
            S.h2o.hm[j, i, k] = (
                S.work1[j0, i0, k] + S.work1[j, i0, k] + S.work1[j2, i0, k] +
                S.work1[j0, i,  k] + S.work1[j, i,  k] + S.work1[j2, i,  k] +
                S.work1[j0, i2, k] + S.work1[j, i2, k] + S.work1[j2, i2, k]
            ) / 9.0f0

            # j = nbneu-1
            j = cfg.nbneu-1
            j0 = j-1
            j2 = j+1
            S.h2o.hm[j, i, k] = (
                S.work1[j0, i0, k] + S.work1[j, i0, k] + S.work1[j2, i0, k] +
                S.work1[j0, i,  k] + 3.0f0 * S.work1[j, i, k] + S.work1[j2, i, k] +
                S.work1[j0, i2, k] + S.work1[j, i2, k] + S.work1[j2, i2, k]
            ) / 11.0f0

            # j = nbneu-2
            j = cfg.nbneu-2
            j0 = j-1
            j2 = j+1
            S.h2o.hm[j, i, k] = (
                S.work1[j0, i0, k] + S.work1[j, i0, k] + S.work1[j2, i0, k] +
                S.work1[j0, i,  k] + 9.0f0 * S.work1[j, i, k] + S.work1[j2, i, k] +
                S.work1[j0, i2, k] + S.work1[j, i2, k] + S.work1[j2, i2, k]
            ) / 17.0f0

            # j = nbneu-3
            j = cfg.nbneu-3
            j0 = j-1
            j2 = j+1
            S.h2o.hm[j, i, k] = (
                S.work1[j0, i0, k] + S.work1[j, i0, k] + S.work1[j2, i0, k] +
                S.work1[j0, i,  k] + 27.0f0 * S.work1[j, i, k] + S.work1[j2, i, k] +
                S.work1[j0, i2, k] + S.work1[j, i2, k] + S.work1[j2, i2, k]
            ) / 35.0f0

            # j = nbneu-4
            j = cfg.nbneu-4
            j0 = j-1
            j2 = j+1
            S.h2o.hm[j, i, k] = (
                S.work1[j0, i0, k] + S.work1[j, i0, k] + S.work1[j2, i0, k] +
                S.work1[j0, i,  k] + 81.0f0 * S.work1[j, i, k] + S.work1[j2, i, k] +
                S.work1[j0, i2, k] + S.work1[j, i2, k] + S.work1[j2, i2, k]
            ) / 89.0f0
        end
    end

    nothing
end


"""
    update_mesopause!(S::MIMASState, cfg::Config)

Recompute the mesopause index `j_mesopause` for each latitude using the current
temperature field `S.uvwtd.tm`. This is called inside the inner timestep loop.
"""
function update_mesopause!(S, cfg)
    nbneu = cfg.nbneu
    kgitneu = cfg.kgitneu
    igitneu = cfg.igitneu

    for j in 1:nbneu
        tmean = zeros(Float32, kgitneu)
        for k in 1:kgitneu
            sum_t = 0.0f0
            for i in 1:igitneu
                sum_t += S.uvwtd.tm[j, i, k]
            end
            tmean[k] = sum_t / igitneu
        end
        # Find level with minimum temperature
        tmin_val, kmeso = findmin(tmean)
        if kmeso > 122
            kmeso = 122
        end
        S.j_mesopause[j] = kmeso
    end
    nothing
end

"""
    particle_birth_control!(S::MIMASState, cfg::Config, ntg)

Fortran equivalent: the two loops (2001 and 2002) after sub_tracertransp.
Resets dust particles that are out of bounds or sporadically repositions them
as a source of new dust. Uses pre‑computed arrays bbreite, bsum, sfactor, invfactor,
and j_mesopause, j_zpress86km, j_zbound_ob.
"""
function particle_birth_control!(S, cfg, ntg)
    # Constants
    bgrenze = 16.0f0   # latitude boundary (yfeld ≤ 16 -> reinitialise)
    nbneu = cfg.nbneu
    igitneu = cfg.igitneu
    kgitneu = cfg.kgitneu
    ntrac = cfg.ntrac

    # Counters (optional, for debugging)
    nbabies = 0
    # nloka1, nloka2, nloka3, ngrup are diagnostic – we can skip or implement later.

    # ------------------------------------------------------------------
    # First loop (2001): every 120th particle (with random offset) if it is dust,
    # reposition it to high latitude and height.
    # ------------------------------------------------------------------
    n_offset = rand(1:120)   # random offset each timestep (like Fortran's int(rand()*120)+1)
    for n in n_offset:120:ntrac
        if S.eis.rfeld[n] == S.eis.rinit[n]   # dust particle
            # xfeld unchanged

            # Latitude selection using the same cumulative distribution as init_staub!
            rando = rand(Float32) * 100.0f0   # 0..100
            ii = 1
            while ii <= 35 && rando > S.bsum[ii]
                ii += 1
            end
            ii = min(ii, 35)
            rando1 = rand(Float32)
            S.eis.yfeld[n] = S.bbreite[ii] + rando1
            if S.eis.yfeld[n] > 53.9999f0
                S.eis.yfeld[n] = 53.999f0
            end

            jj = floor(Int, S.eis.yfeld[n])
            jj = clamp(jj, 1, nbneu)

            rando = rand(Float32)
            if rando > 0.8f0
                # Upper tail (20%)
                rando2 = (rando - 0.8f0) * 5.0f0
                kk1 = S.j_mesopause[jj] - (5 + S.j_zpress86km[jj])
                kk2 = 15
                if kk1 >= 1
                    kk2 += kk1
                end
                S.eis.zfeld[n] = Float32(S.j_zpress86km[jj]) - 3.0f0 + Float32(kk2) * rando2
            else
                # Main 80%
                rando2 = rando * (1.0f0 / 0.8f0)
                S.eis.zfeld[n] = Float32(S.j_zpress86km[jj]) - 5.0f0 + 10.0f0 * rando2
            end

            # Clamp z
            if S.eis.zfeld[n] > kgitneu
                S.eis.zfeld[n] = kgitneu - eps
            elseif S.eis.zfeld[n] < 1.0f0
                S.eis.zfeld[n] = 1.0f0 + eps
            end
        end
    end

    # ------------------------------------------------------------------
    # Second loop (2002): all particles – if dust and out of bounds, reset.
    # ------------------------------------------------------------------
    for n in 1:ntrac
        if S.eis.rfeld[n] == S.eis.rinit[n]   # dust particle
            k = floor(Int, S.eis.zfeld[n])
            j = floor(Int, S.eis.yfeld[n])
            if j > nbneu
                # Should not happen, but for safety
                j = nbneu
            elseif j < 1
                j = 1
            end
            if S.eis.yfeld[n] <= bgrenze || k <= 62 || k > S.j_zbound_ob[j]
                nbabies += 1

                # Reset radius to initial
                S.eis.rfeld[n] = S.eis.rinit[n]

                # Longitude: random in [1, igitneu]
                rando = rand(Float32)
                S.eis.xfeld[n] = 1.0f0 + rando * (igitneu - 1.0f0)
                if S.eis.xfeld[n] > igitneu
                    S.eis.xfeld[n] = igitneu - eps
                elseif S.eis.xfeld[n] < 1.0f0
                    S.eis.xfeld[n] = 1.0f0 + eps
                end

                # Latitude: same distribution as first loop
                rando = rand(Float32) * 100.0f0
                ii = 1
                while ii <= 35 && rando > S.bsum[ii]
                    ii += 1
                end
                ii = min(ii, 35)
                rando1 = rand(Float32)
                S.eis.yfeld[n] = S.bbreite[ii] + rando1
                if S.eis.yfeld[n] > 53.9999f0
                    S.eis.yfeld[n] = 53.999f0
                end

                jj = floor(Int, S.eis.yfeld[n])
                jj = clamp(jj, 1, nbneu)

                rando = rand(Float32)
                if rando > 0.8f0
                    rando2 = (rando - 0.8f0) * 5.0f0
                    kk1 = S.j_mesopause[jj] - (5 + S.j_zpress86km[jj])
                    kk2 = 15
                    if kk1 >= 1
                        kk2 += kk1
                    end
                    S.eis.zfeld[n] = Float32(S.j_zpress86km[jj]) - 3.0f0 + Float32(kk2) * rando2
                else
                    rando2 = rando * (1.0f0 / 0.8f0)
                    S.eis.zfeld[n] = Float32(S.j_zpress86km[jj]) - 5.0f0 + 10.0f0 * rando2
                end

                if S.eis.zfeld[n] > kgitneu
                    S.eis.zfeld[n] = kgitneu - eps
                elseif S.eis.zfeld[n] < 1.0f0
                    S.eis.zfeld[n] = 1.0f0 + eps
                end
            end
        end
    end

    # Optional: print number of resets occasionally (for debugging)
    # if nbabies > 0
    #     println("        nbabies = ", nbabies)
    # end
    nothing
end

# ----------------------------------------------------------------------
# 1D advection routine (sub_advec1d) – corrected version
# ----------------------------------------------------------------------
"""
    sub_advec1d!(ix, iperiod, dxx, u, den0, den1, dd0, q0, qn, cfg)

Fortran equivalent: `sub_advec1d`.
Performs 1D advection with monotonic flux limiting.
- ix: number of interior cells
- iperiod: 1 if periodic (zonal), 0 otherwise (meridional/vertical)
- dxx: grid spacing (scalar)
- u: velocity at cell faces (0:ix) – stored in a vector of length ix+1,
      where u[1] corresponds to u(0), u[ix+1] to u(ix)
- den0: density at beginning of step (1:ix)
- den1: density at end of step (1:ix)
- dd0: density at cell faces (0:ix) – length ix+1, dd0[1] = dd0(0)
- q0: input mixing ratio with halo (0:ix+1) – length ix+2, q0[1] = q0(0)
- qn: output mixing ratio (1:ix) – length ix
- cfg: configuration (provides dttrans)
"""
function sub_advec1d!(ix, iperiod, dxx, u, den0, den1, dd0, q0, qn, cfg)
    # flux at faces 0..ix (size ix+1)
    flux = zeros(Float32, ix+1)
    vcmax = zeros(Float32, ix)
    vcmin = zeros(Float32, ix)
    imxmn = falses(ix+2)   # indices 1..ix+2 correspond to Fortran 0..ix+1

    zr0 = 0.0f0

    # Identify local maxima/minima for interior cells (i=1..ix in Fortran)
    for i in 1:ix
        # Fortran index i corresponds to Julia q0 index i+1
        ii = i + 1
        imxmn[ii] = (q0[ii-1] - q0[ii]) * (q0[ii] - q0[ii+1]) <= zr0

        ck1 = q0[ii]
        ck2 = q0[ii]
        if u[i+1] < zr0   # u(i) < 0? Wait: u(i) is u[i+1] in Julia
            # Actually u(i) is u[i+1] because u[1]=u(0). So u(i) < 0 => u[i+1] < 0
            ck1 = q0[ii+1]
        end
        if u[i] > zr0      # u(i-1) > 0? u(i-1) is u[i] (since u[1]=u(0))
            ck2 = q0[ii-1]
        end
        vcmax[i] = max(q0[ii], ck1, ck2)
        vcmin[i] = min(q0[ii], ck1, ck2)
    end

    # Upwind flux at left boundary if inflow (u(0) >= 0)
    if u[1] >= zr0
        flux[1] = q0[1] * u[1] * cfg.dttrans * dd0[1]
    end

    # Forward sweep (cells with u(i) >= 0)
    for i in 1:ix
        if u[i+1] >= zr0   # u(i) >= 0
            if u[i] < zr0   # u(i-1) < 0 (outflow-only cell)
                flux[i+1] = q0[i+1] * u[i+1] * cfg.dttrans * dd0[i+1]
            else
                x1 = u[i+1] * cfg.dttrans / dxx
                x1n = (1.0f0 - x1) * (q0[i+2] - q0[i]) / 4.0f0
                cf = q0[i+1] + x1n

                if imxmn[i]   # i-1 (since imxmn stored for index i+1? Wait careful: imxmn[ii] corresponds to Fortran i, so for current i, we need imxmn for i-1 and i+1. imxmn indices: imxmn[i+1] is for cell i (Fortran). So imxmn[i] is for cell i-1, imxmn[i+2] is for cell i+1.
                    # Let's define: imxmn[ii] corresponds to Fortran cell i. So for current i, imxmn[i] is cell i-1, imxmn[i+2] is cell i+1.
                    cf = q0[i+1] + max(1.5f0, 1.2f0 + 0.6f0 * x1) * x1n
                end
                if imxmn[i+2]
                    cf = q0[i+1] + (1.75f0 - 0.45f0 * x1) * x1n
                end

                cf1 = clamp(cf, min(q0[i+1], q0[i+2]), max(q0[i+1], q0[i+2]))

                qn[i] = max(vcmin[i], min(vcmax[i],
                    (q0[i+1] * den0[i] - x1 * cf1 * dd0[i+1] + flux[i] / dxx) / den1[i]
                ))

                flux[i+1] = dxx * (q0[i+1] * den0[i] - qn[i] * den1[i]) + flux[i]
            end
        end
    end

    # Periodic boundary adjustment (forward) – if iperiod == 1
    if iperiod == 1
        if u[ix] >= zr0 && u[ix+1] >= zr0   # u(idim-1) and u(idim)
            qn[1] = (q0[2] * den0[1] - flux[2] / dxx + flux[ix+1] / dxx) / den1[1]
        end
    end

    # Downwind sweep for u < 0
    if u[ix+1] < zr0   # u(ix) < 0
        flux[ix+1] = q0[ix+2] * u[ix+1] * cfg.dttrans * dd0[ix+1]
    end

    for i in ix:-1:1
        if u[i] >= zr0   # u(i-1) >= 0
            if u[i+1] < zr0   # u(i) < 0 (inflow-only cell)
                qn[i] = max(vcmin[i], min(vcmax[i],
                    (q0[i+1] * den0[i] - flux[i+1] / dxx + flux[i] / dxx) / den1[i]
                ))
            end
        else
            x1 = cfg.dttrans * abs(u[i]) / dxx   # u(i-1) is u[i]
            x1n = (1.0f0 - x1) * (q0[i] - q0[i+2]) / 4.0f0
            cf = q0[i+1] + x1n
            if imxmn[i+2]   # i+1
                cf = q0[i+1] + max(1.5f0, 1.2f0 + 0.6f0 * x1) * x1n
            end
            if imxmn[i]     # i-1
                cf = q0[i+1] + (1.75f0 - 0.45f0 * x1) * x1n
            end
            cf1 = clamp(cf, min(q0[i+1], q0[i]), max(q0[i+1], q0[i]))
            if u[i+1] >= zr0   # u(i) >= 0
                cf1 = q0[i+1]   # outflow-only cell
            end
            qn[i] = max(vcmin[i], min(vcmax[i],
                (q0[i+1] * den0[i] - flux[i+1] / dxx - x1 * cf1 * dd0[i]) / den1[i]
            ))
            flux[i] = dxx * (qn[i] * den1[i] - q0[i+1] * den0[i]) + flux[i+1]
        end
    end

    # Periodic boundary adjustment (backward)
    if iperiod == 1
        if u[2] < zr0 && u[ix+1] < zr0   # u(1) and u(idim)
            qn[ix] = (q0[ix+1] * den0[ix] - flux[1] / dxx + flux[ix] / dxx) / den1[ix]
        end
    end

    nothing
end


# ----------------------------------------------------------------------
# 3D advection (sub_transwalcek) – corrected version using constants dy, dz
# ----------------------------------------------------------------------
function sub_transwalcek!(S, cfg)
    d0 = 1.0f0

    # Zonal advection
    iperiod = 1
    for k in 22:cfg.kgitneu-2
        for j in 2:cfg.nbneu
            deltax = S.consts.dx[j]
            q0i = zeros(Float32, cfg.igitneu+2)
            ui = zeros(Float32, cfg.igitneu+1)
            den0i = zeros(Float32, cfg.igitneu)
            den1i = zeros(Float32, cfg.igitneu)
            dd0i = zeros(Float32, cfg.igitneu+1)

            for i in 1:cfg.igitneu
                q0i[i+1] = S.h2o.hm[j, i, k]
            end
            q0i[1] = q0i[cfg.igitneu+1]
            q0i[cfg.igitneu+2] = q0i[2]

            for i in 1:cfg.igitneu
                ui[i] = S.uvwtd.um[j, i, k]
            end
            ui[cfg.igitneu+1] = ui[1]

            for i in 1:cfg.igitneu
                den0i[i] = d0
                den1i[i] = den0i[i] - cfg.dttrans / deltax * (d0 * ui[i+1] - d0 * ui[i])
            end
            dd0i .= d0

            qni = zeros(Float32, cfg.igitneu)
            sub_advec1d!(cfg.igitneu, iperiod, deltax, ui, den0i, den1i, dd0i, q0i, qni, cfg)

            for i in 1:cfg.igitneu
                S.h2o.hm[j, i, k] = qni[i]
            end
        end
    end

    # Meridional advection
    iperiod = 0
    for k in 22:cfg.kgitneu-2
        for i in 1:cfg.igitneu
            deltax = dy   # use constant dy from constants.jl
            q0j = zeros(Float32, cfg.nbneu+2)
            uj = zeros(Float32, cfg.nbneu+1)
            den0j = zeros(Float32, cfg.nbneu)
            den1j = zeros(Float32, cfg.nbneu)
            dd0j = zeros(Float32, cfg.nbneu+1)

            for j in 1:cfg.nbneu
                q0j[j+1] = S.h2o.hm[j, i, k]
            end
            q0j[1] = q0j[2]
            q0j[cfg.nbneu+2] = q0j[cfg.nbneu+1]

            for j in 1:cfg.nbneu
                uj[j] = S.uvwtd.vm[j, i, k]
            end
            uj[cfg.nbneu+1] = S.uvwtd.vm[cfg.nbneu+1, i, k]

            ip1 = i == cfg.igitneu ? 1 : i+1
            im1 = i == 1 ? cfg.igitneu : i-1
            for j in 1:cfg.nbneu
                den0j[j] = d0 - cfg.dttrans / S.consts.dx[j] *
                           (d0 * S.uvwtd.um[j, ip1, k] - d0 * S.uvwtd.um[j, im1, k])
                den1j[j] = den0j[j] - cfg.dttrans / deltax *
                           (d0 * uj[j+1] - d0 * uj[j])
            end
            dd0j .= d0

            qnj = zeros(Float32, cfg.nbneu)
            sub_advec1d!(cfg.nbneu, iperiod, deltax, uj, den0j, den1j, dd0j, q0j, qnj, cfg)

            for j in 2:cfg.nbneu
                S.h2o.hm[j, i, k] = qnj[j]
            end
        end
    end

    # Vertical advection
    iperiod = 0
    klevels = cfg.kgitneu - cfg.lunten - 1
    klevels1 = klevels + 1
    for i in 1:cfg.igitneu
        for j in 2:cfg.nbneu
            deltax = dz   # use constant dz from constants.jl
            q0k = zeros(Float32, klevels1+1)
            uk = zeros(Float32, klevels+1)
            den0k = zeros(Float32, klevels)
            den1k = zeros(Float32, klevels)
            dd0k = zeros(Float32, klevels+1)

            for kk in 1:klevels1
                k = kk + cfg.lunten - 1
                q0k[kk+1] = S.h2o.hm[j, i, k]
            end
            q0k[1] = S.h2o.hm[j, i, cfg.lunten-1]

            for kk in 0:klevels
                uk[kk+1] = S.uvwtd.wm[j, i, kk + cfg.lunten]
            end

            ip1 = i == cfg.igitneu ? 1 : i+1
            im1 = i == 1 ? cfg.igitneu : i-1
            for kk in 1:klevels
                k = kk + cfg.lunten - 1
                den0k[kk] = d0 -
                    cfg.dttrans / S.consts.dx[j] * (d0 * S.uvwtd.um[j, ip1, k] - d0 * S.uvwtd.um[j, im1, k]) -
                    cfg.dttrans / dy * (d0 * S.uvwtd.vm[j+1, i, k] - d0 * S.uvwtd.vm[j, i, k])
                den1k[kk] = den0k[kk] - cfg.dttrans / deltax *
                            (d0 * uk[kk+1] - d0 * uk[kk])
            end
            dd0k .= d0

            qnk = zeros(Float32, klevels)
            sub_advec1d!(klevels, iperiod, deltax, uk, den0k, den1k, dd0k, q0k, qnk, cfg)

            for kk in 1:klevels
                k = kk + cfg.lunten - 1
                S.h2o.hm[j, i, k] = qnk[kk]
            end
        end
    end

    nothing
end



# ----------------------------------------------------------------------
# Tridiagonal solver (tridiag)
# ----------------------------------------------------------------------
"""
    tridiag!(a, b, c, r, u, n)

Fortran equivalent: `tridiag`.
Solves a tridiagonal system. a, b, c are the sub‑, main, and super‑diagonals.
r is the right‑hand side, u is the solution. n is the system size.
Uses the Thomas algorithm.
"""
function tridiag!(a, b, c, r, u, n)
    gam = zeros(Float32, n)
    if b[1] == 0.0f0
        error("tridiag: zero pivot")
    end
    bet = b[1]
    u[1] = r[1] / bet

    for j in 2:n
        gam[j] = c[j-1] / bet
        bet = b[j] - a[j] * gam[j]
        if bet == 0.0f0
            error("tridiag: zero pivot")
        end
        u[j] = (r[j] - a[j] * u[j-1]) / bet
    end

    for j in n-1:-1:1
        u[j] -= gam[j+1] * u[j+1]
    end
    nothing
end
# ----------------------------------------------------------------------
# Vertical diffusion of H2O (sub_diffu_h2o)
# ----------------------------------------------------------------------
function sub_diffu_h2o!(S, cfg)
    dttrans_2 = cfg.dttrans_2
    for i in 1:cfg.igitneu
        for j in 3:cfg.nbneu   # start from j=3 as in Fortran
            a = zeros(Float32, cfg.kgitneu)
            b = zeros(Float32, cfg.kgitneu)
            c = zeros(Float32, cfg.kgitneu)
            r = zeros(Float32, cfg.kgitneu)
            u = zeros(Float32, cfg.kgitneu)

            for k in 1:cfg.kgitneu
                alpha = S.consts.turbkzz[k] * dttrans_2 / (dz * dz)
                a[k] = -alpha
                b[k] = 1.0f0 + 2.0f0 * alpha
                c[k] = -alpha
                r[k] = S.h2o.hm[j, i, k]
            end
            a[1] = 0.0f0
            c[cfg.kgitneu] = 0.0f0

            tridiag!(a, b, c, r, u, cfg.kgitneu)

            for k in cfg.lunten:cfg.kgitneu-2
                S.h2o.hm[j, i, k] = u[k]
            end
        end
    end
    nothing
end