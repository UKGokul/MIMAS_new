# src/init.jl
# Initialization of the MIMAS model state.
# All functions follow the Fortran logic and print progress for transparency.
# This file contains routines that fill the model state with geometry,
# derived indices, lookup tables, and other constant data needed before
# the main time loop begins.

include("constants.jl")
include("config.jl")
include("types.jl")

# ----------------------------------------------------------------------
# Step E1: Allocate all arrays with zeros (safe initialisation)
# ----------------------------------------------------------------------
"""
    allocate_state(cfg::Config) -> MIMASState

Allocate all arrays of the model state with zeros.
Dimensions are taken from the configuration.
This is the first step after defining the configuration – it creates
a fresh state object where all fields are zeroed.
"""
allocate_state(cfg::Config) = MIMASState(cfg)


# ----------------------------------------------------------------------
# Step E2: Fill constant geometry arrays (dx, zgeo) and read fixed files
# ----------------------------------------------------------------------
"""
    init_constants!(S::MIMASState, cfg::Config)

Fill geometry arrays `dx` (zonal grid spacing) and `zgeo` (geometric height).
`dx` is computed from latitude and Earth radius, `zgeo` is a linear height
array from 77.8 km upward in 0.1 km steps.
This function is called immediately after state allocation.
"""
function init_constants!(S, cfg)
    println("  Computing geometry arrays...")

    # dx: zonal grid spacing (depends on latitude)
    # In Fortran: dx(j) = 2π/igitneu * rearth * cos(lat) with lat in degrees.
    for j in 1:cfg.nbneu
        lat = 37.5 + (j-1)   # latitude in degrees (as in Fortran)
        S.consts.dx[j] = Float32(2π / cfg.igitneu * rearth * cosd(lat))
    end

    # zgeo: geometric height (km) from 77.8 km upward, step 0.1 km
    for k in 1:cfg.kgitneu
        S.consts.zgeo[k] = 77.8f0 + (k-1) * 0.1f0
    end

    # Print first value and range for verification
    println("    dx[1] = ", S.consts.dx[1], " m, dx[end] = ", S.consts.dx[end], " m")
    println("    zgeo[1] = ", S.consts.zgeo[1], " km, zgeo[end] = ", S.consts.zgeo[end], " km")
end


# ----------------------------------------------------------------------
# Step E3: Compute mesopause and 86 km pressure level indices
# ----------------------------------------------------------------------
"""
    compute_mesopause_and_86km!(S::MIMASState, cfg::Config)

For each latitude, compute the mean temperature and pressure altitude profiles,
find the mesopause level (minimum temperature, capped at 122) and the level
where pressure altitude crosses 89.1 km (for 86 km geometric height reference).
Stores results in `S.j_mesopause` and `S.j_zpress86km`.
These indices are used later for particle initialisation and birth control.
"""
function compute_mesopause_and_86km!(S, cfg)
    nbneu = cfg.nbneu
    kgitneu = cfg.kgitneu
    igitneu = cfg.igitneu

    for j in 1:nbneu
        tmean = zeros(Float32, kgitneu)
        zpressmean = zeros(Float32, kgitneu)
        for k in 1:kgitneu
            sum_t = 0.0f0
            sum_z = 0.0f0
            for i in 1:igitneu
                sum_t += S.uvwtd48.tm1[j, i, k]
                sum_z += S.uvwtd48.zm1[j, i, k]
            end
            tmean[k] = sum_t / igitneu
            zpressmean[k] = sum_z / igitneu
        end

        # Find level with minimum temperature (mesopause)
        tmin_val, kmeso = findmin(tmean)
        if kmeso > 122
            kmeso = 122
        end
        S.j_mesopause[j] = kmeso

        # Find level where zpressmean crosses 89.1 km
        k86 = 0
        for k in 1:kgitneu-1
            if zpressmean[k] <= 89.1f0 && zpressmean[k+1] >= 89.1f0
                k86 = k
                break
            end
        end
        S.j_zpress86km[j] = k86
    end

    # Print a sample for first and last latitude to verify
    println("  Computed mesopause and 86 km indices.")
    println("    j_mesopause[1] = ", S.j_mesopause[1], ", j_mesopause[end] = ", S.j_mesopause[end])
    println("    j_zpress86km[1] = ", S.j_zpress86km[1], ", j_zpress86km[end] = ", S.j_zpress86km[end])
end


# ----------------------------------------------------------------------
# Step E4: Compute all lookup tables (microphysics, sedimentation, etc.)
# ----------------------------------------------------------------------
"""
    init_tables!(S::MIMASState, cfg::Config)

Compute all pre‑computed tables needed during the time loop:
- horiwicht_tab   : horizontal weighting factor for particle‑water coupling.
- vpsatt_tab      : saturation vapour pressure over ice (Murphy & Koop 2005).
- sqrt1t_tab      : factor for ice growth rate.
- sqrtsedi_tab    : factor for sedimentation velocity.
- faktkelvin_tab  : Kelvin (curvature) effect on saturation vapour pressure.
- wrichtig_tab    : vertical random walk amplitude (from turbulent diffusion).
- volufak_tab     : vertical volume scaling factor.
- worka_tab       : particle temperature adjustment factor (radius‑ & height‑dependent).
All tables are stored in `S.tab`.
"""
function init_tables!(S, cfg)
    println("  Computing lookup tables...")

    # Extract constants
    TAB_TMAX = 71
    TAB_RMAX = 300
    TAB_RPTP = 200
    nbneu = cfg.nbneu
    kgitneu = cfg.kgitneu
    igitneu = cfg.igitneu

    # ------------------------------------------------------------------
    # horiwicht_tab (horizontal weighting factor)
    # ------------------------------------------------------------------
    # Formula: horiwicht_tab(j) = 1e12 * igitneu / (cos(lat) * 60)
    # This factor is used in the microphysics to couple particle water exchange.
    for j in 1:nbneu
        lat = Float32(37.5) + Float32(j - 1)
        S.tab.horiwicht_tab[j] =
            Float32(1e12) * igitneu / (cosd(lat) * Float32(60.0))
    end
    println("    horiwicht_tab[1] = ", S.tab.horiwicht_tab[1])

    # ------------------------------------------------------------------
    # Tables that depend on temperature (100..170 K)
    # ------------------------------------------------------------------
    for itemp in 1:TAB_TMAX
        tp = Float32(99.0) + Float32(itemp)   # temperature in K

        # Saturation vapour pressure over ice (Murphy & Koop 2005)
        # Used in microphysics to compute supersaturation.
        S.tab.vpsatt_tab[itemp] = exp(
            Float32(9.550426)
            - Float32(5723.265) / tp
            + Float32(3.53068) * log(tp)
            - Float32(0.00728332) * tp
        )

        # sqrt(1/T) factor for growth rate (drdt)
        S.tab.sqrt1t_tab[itemp] =
            Float32(1e9) * Float32(0.83) *
            sqrt(xmh2o / (Float32(2) * Float32(π) * xkbolz * tp)) /
            rhoice

        # Sedimentation factor (wfall = rfeld * sqrtsedi_tab / dm)
        S.tab.sqrtsedi_tab[itemp] =
            Float32(1e-9) * rhoice * gmeso *
            sqrt(Float32(π) * xkbolz * tp / xmair) /
            (sqrt(Float32(2.0)) * Float32(2.0))

        # Kelvin factor for each radius bin (0.1–30 nm)
        # This accounts for curvature effect on saturation vapour pressure.
        for iradi in 1:TAB_RMAX
            rwert0 = Float32(iradi) * Float32(0.1)   # radius in nm
            sigma = (Float32(0.141) - Float32(1.5e-4) * tp) /
                    (Float32(1.0) + Float32(2.0) * Float32(0.15) / rwert0)
            S.tab.faktkelvin_tab[iradi, itemp] = exp(
                Float32(1e9) * Float32(2.0) * xmh2o * sigma /
                (rhoice * xkbolz * tp * rwert0)
            )
        end
    end

    # Adjust Kelvin factor for continuity: subtract value at 30 nm
    # This ensures a smooth transition between growth and sublimation.
    for itemp in 1:TAB_TMAX
        w30 = S.tab.faktkelvin_tab[TAB_RMAX, itemp] - Float32(1.0)
        for iradi in 1:TAB_RMAX
            S.tab.faktkelvin_tab[iradi, itemp] -= w30
        end
    end

    # Print a few sample values to verify
    println("    vpsatt_tab[1] (T=100K) = ", S.tab.vpsatt_tab[1])
    println("    sqrt1t_tab[1] = ", S.tab.sqrt1t_tab[1])
    println("    sqrtsedi_tab[1] = ", S.tab.sqrtsedi_tab[1])
    println("    faktkelvin_tab[1,1] = ", S.tab.faktkelvin_tab[1,1])

    # ------------------------------------------------------------------
    # Vertical turbulence and volume factors (depend on zgeo)
    # ------------------------------------------------------------------
    # wrichtig_tab is used for random walk amplitude in vertical particle transport.
    # volufak_tab converts particle counts to number density.
    for k in 1:kgitneu
        S.tab.wrichtig_tab[k] = sqrt(S.consts.turbkzz[k] / Float32(7.5))
        S.tab.volufak_tab[k] =
            Float32(1.0) / exp((Float32(88.0) - S.consts.zgeo[k]) / Float32(4.0))
    end
    println("    wrichtig_tab[1] = ", S.tab.wrichtig_tab[1])
    println("    volufak_tab[1] = ", S.tab.volufak_tab[1])

    # ------------------------------------------------------------------
    # Particle temperature factor (worka_tab)
    # ------------------------------------------------------------------
    # This factor modifies the particle temperature based on radius and height.
    # For radii ≤10 nm: scale = 0; for radii ≥30 nm: scale = 1;
    # linear in between. Then multiplied by an exponential height factor.
    for k in 1:kgitneu
        for iradi in 1:TAB_RPTP
            scale = iradi <= 10 ? Float32(0.0) :
                    iradi >= 30 ? Float32(1.0) :
                    Float32(iradi - 10) / Float32(20.0)

            S.tab.worka_tab[iradi, k] =
                scale *
                Float32(iradi) *
                Float32(0.01) *
                beta05 *
                t82km *
                exp((S.consts.zgeo[k] - Float32(82.0)) / Float32(4.5))
        end
    end
    println("    worka_tab[10,1] = ", S.tab.worka_tab[10,1])
    println("    worka_tab[20,1] = ", S.tab.worka_tab[20,1])
    println("    worka_tab[30,1] = ", S.tab.worka_tab[30,1])

    println("  init_tables! done.")
end

"""
    compute_j_zbound_ob!(S::MIMASState, cfg::Config)

Compute the upper boundary index `j_zbound_ob` for each latitude, used in particle birth control.
This follows the Fortran loop that sets values for j=1..10 and spreads them over latitudes.
The array is constant for a given grid.
"""
function compute_j_zbound_ob!(S, cfg)
    nbneu = cfg.nbneu
    S.j_zbound_ob .= 0
    # Set base values for the top latitudes (j around 44..53 etc.)
    for jj in 1:10
        idx = 43 + jj
        S.j_zbound_ob[idx] = 141        # +800 m
        idx = 33 + jj
        S.j_zbound_ob[idx] = 130 + jj
        idx = 23 + jj
        S.j_zbound_ob[idx] = 120 + jj
        idx = 13 + jj
        S.j_zbound_ob[idx] = 110 + jj
        idx = 3 + jj
        S.j_zbound_ob[idx] = 100 + jj
    end
    # Ensure indices are within 1..nbneu (they are, given nbneu=53)
    println("    j_zbound_ob computed (sample: first 5 = ", S.j_zbound_ob[1:5], ")")
    nothing
end