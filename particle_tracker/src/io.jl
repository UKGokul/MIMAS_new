# src/io.jl
# File I/O routines for MIMAS
# This file contains all functions that read from or write to external files.
# It handles both formatted (text) and unformatted (binary) Fortran-style files,
# with proper endian conversion (big‑endian as used in the original model).
# The main subroutines included here are:
#   - read_fortran_record: low‑level reader for Fortran unformatted records.
#   - init_global!: reads multiple global initialisation files (turbulence, backscatter,
#                   Earth ellipse, Lyman‑alpha, UT1, CH4).
#   - read_h2o_init!: reads the H2O pressure profile file.
#   - read_dynamics!: reads a full dynamics snapshot (including mean‑file corrections)
#                     into the *1 arrays.
#   - write_fortran_record: low‑level writer for Fortran unformatted records.
#   - write_restart!: writes the restart file with particle data.

using Sockets: ntoh   # for endian conversion (network to host)
using NetCDF
# ----------------------------------------------------------------------
# Low‑level helper: read a Fortran unformatted record
# ----------------------------------------------------------------------
"""
    read_fortran_record(io::IO, T::Type, dims...)

Read one Fortran unformatted record from `io`. Returns an array of type `T`
with dimensions `dims`. Handles the 4‑byte header/trailer and big‑endian conversion.
"""
function read_fortran_record(io::IO, T::Type, dims...)
    # Each Fortran unformatted record starts with a 4‑byte integer giving the
    # number of bytes in the record (big‑endian). We read it, convert to native order.
    header = read(io, Int32)
    header = ntoh(header)   # convert from big‑endian to native

    # Compute how many bytes we expect based on requested dimensions.
    expected_bytes = sizeof(T) * prod(dims)
    if header != expected_bytes
        error("Record size mismatch: expected $expected_bytes bytes, got $header")
    end

    # Allocate an uninitialised array of the correct shape and read the data.
    data = Array{T}(undef, dims)
    read!(io, data)

    # If the machine is little‑endian (most common), the bytes read are in
    # big‑endian order and must be swapped. bswap works for Float32 as well.
    if T <: Number
        for i in eachindex(data)
            data[i] = bswap(data[i])
        end
    end

    # Read the trailer (another 4‑byte integer, should equal the header).
    trailer = read(io, Int32)
    trailer = ntoh(trailer)
    if trailer != header
        error("Trailer does not match header: $trailer vs $header")
    end

    return data
end

# ----------------------------------------------------------------------
# Global initialisation (sub_init_global)
# ----------------------------------------------------------------------
"""
    init_global!(S::MIMASState, cfg::Config, year::Int) -> (Vector{Float32}, Float32, Float32)

Fortran equivalent: `sub_init_global` (called early in the main program).
Reads all global initialisation files that are independent of the particular
simulation date, except for the Lyman‑α file which is year‑dependent.
Fills the corresponding arrays in `S` and returns:
- xlyman_obs: daily Lyman‑α flux values (length 366, only first max_doy used)
- scale_ch4: CH4 scaling factor for H2O (fraction, multiplied later)
- unix_time: reference Unix time (seconds since 1970-01-01) for 10 May of the given year
The arrays filled in S.consts and S.h2o are:
    turbkzz      (vertical diffusion coefficient)
    backradius, crosssec*, ext*  (optical properties for particles)
    ellipsefak, oblinoa, true_solar_time (Earth orbit parameters)
"""
function init_global!(S::MIMASState, cfg::Config, year::Int)
    println("  Initialising global data from files...")

    # Helper for leap year (as in Fortran, except 1900 is handled specially)
    function days_in_year(y)
        if y == 1900
            return 365
        else
            return (y % 4 == 0) ? 366 : 365
        end
    end
    max_doy = days_in_year(year)

    # ------------------------------------------------------------------
    # 1. Turbulence profile (CTURBKZZ) – binary, one record containing two 286‑element arrays
    # ------------------------------------------------------------------
    println("    Reading turbulence profile: $(cfg.cturbkzz)")
    try
        open(cfg.cturbkzz, "r") do io
            # Read one record containing 572 floats (286 * 2)
            data = read_fortran_record(io, Float32, 572)
            turbkzz2d = data[1:286]
            wturb2d   = data[287:572]   # not used further, but read for completeness

            # Interpolate (every second element) and scale by 0.25
            k2 = 0
            for k in 2:2:286
                k2 += 1
                if k2 <= cfg.kgitneu
                    S.consts.turbkzz[k2] = turbkzz2d[k] * 0.25f0
                end
            end

            # Extrapolate top levels (k=144 to 163 in Fortran indexing)
            for kk in 1:20
                k = 143 + kk
                if k <= cfg.kgitneu
                    faktkzz = 1.0f0 + Float32(kk) / 20.0f0
                    S.consts.turbkzz[k] = S.consts.turbkzz[143] / faktkzz
                end
            end

            # Set lower levels to constant (Fortran modification)
            k2_override = 23 + 50   # as in Fortran comment
            for k in 1:min(k2_override, cfg.kgitneu)
                if k2_override+1 <= cfg.kgitneu
                    S.consts.turbkzz[k] = S.consts.turbkzz[k2_override+1]
                end
            end
        end
        # Print first and last few values of turbkzz to verify
        println("      turbkzz first 5: ", S.consts.turbkzz[1:5])
        println("      turbkzz last 5:  ", S.consts.turbkzz[end-4:end])
    catch e
        println("      Warning: could not read turbulence file, using dummy values (all 10.0)")
        fill!(S.consts.turbkzz, 10.0f0)
    end

    # ------------------------------------------------------------------
    # 2. Backscatter data (CBACKSCATTER532) – binary, one record containing 11 * 1500 floats
    # ------------------------------------------------------------------
    println("    Reading backscatter data: $(cfg.cbackscatter532)")
    try
        open(cfg.cbackscatter532, "r") do io
            data = read_fortran_record(io, Float32, 11 * 1500)  # 16500 floats
            offset = 0
            S.consts.backradius   .= data[offset+1:offset+1500]; offset += 1500
            S.consts.crosssec532  .= data[offset+1:offset+1500]; offset += 1500
            S.consts.ext532       .= data[offset+1:offset+1500]; offset += 1500
            S.consts.crosssec126  .= data[offset+1:offset+1500]; offset += 1500
            S.consts.ext126       .= data[offset+1:offset+1500]; offset += 1500
            S.consts.crosssec200  .= data[offset+1:offset+1500]; offset += 1500
            S.consts.ext200       .= data[offset+1:offset+1500]; offset += 1500
            S.consts.crosssec1000 .= data[offset+1:offset+1500]; offset += 1500
            S.consts.ext1000      .= data[offset+1:offset+1500]; offset += 1500
            S.consts.crosssec3000 .= data[offset+1:offset+1500]; offset += 1500
            S.consts.ext3000      .= data[offset+1:offset+1500]
        end
        # Print first and last few of backradius and crosssec532 to verify
        println("      backradius first 5: ", S.consts.backradius[1:5])
        println("      backradius last 5:  ", S.consts.backradius[end-4:end])
        println("      crosssec532 first 5: ", S.consts.crosssec532[1:5])
        println("      crosssec532 last 5:  ", S.consts.crosssec532[end-4:end])
    catch e
        println("      Warning: could not read backscatter file, using dummy zeros")
        fill!(S.consts.backradius, 0.0f0)
        fill!(S.consts.crosssec532, 0.0f0)
        fill!(S.consts.ext532, 0.0f0)
        fill!(S.consts.crosssec126, 0.0f0)
        fill!(S.consts.ext126, 0.0f0)
        fill!(S.consts.crosssec200, 0.0f0)
        fill!(S.consts.ext200, 0.0f0)
        fill!(S.consts.crosssec1000, 0.0f0)
        fill!(S.consts.ext1000, 0.0f0)
        fill!(S.consts.crosssec3000, 0.0f0)
        fill!(S.consts.ext3000, 0.0f0)
    end

    # ------------------------------------------------------------------
    # 3. Earth ellipse data (CELLIPSE_EARTH) – 365 separate records, each with 4 floats
    # ------------------------------------------------------------------
    println("    Reading Earth ellipse data: $(cfg.cellipse_earth)")
    try
        open(cfg.cellipse_earth, "r") do io
            for i in 1:365
                data = read_fortran_record(io, Float32, 4)   # read 4 floats per record
                S.h2o.ellipsefak[i]       = data[2]
                S.h2o.oblinoa[i]          = data[3]
                S.h2o.true_solar_time[i]  = data[4]
            end
        end
        println("      ellipsefak first 5: ", S.h2o.ellipsefak[1:5])
        println("      ellipsefak last 5:  ", S.h2o.ellipsefak[end-4:end])
        println("      oblinoa first 5: ", S.h2o.oblinoa[1:5])
        println("      oblinoa last 5:  ", S.h2o.oblinoa[end-4:end])
        println("      true_solar_time first 5: ", S.h2o.true_solar_time[1:5])
        println("      true_solar_time last 5:  ", S.h2o.true_solar_time[end-4:end])
    catch e
        println("      Warning: could not read Earth ellipse file, using dummy values")
        fill!(S.h2o.ellipsefak, 1.0f0)
        fill!(S.h2o.oblinoa, 23.5f0)
        fill!(S.h2o.true_solar_time, 0.0f0)
    end

    # ------------------------------------------------------------------
    # 4. Lyman‑alpha data – formatted file per year
    # ------------------------------------------------------------------
    println("    Reading Lyman‑alpha data for year $year")
    lyman_base = hasproperty(cfg, :lyman_pfad) ? cfg.lyman_pfad : "/home/gu01/HPC/klymanalpha/"
    lyman_file = joinpath(lyman_base, "$(lpad(year,4,'0')).lyman_alpha.txt")
    xlyman_obs = zeros(Float32, 366)
    try
        open(lyman_file, "r") do io
            for i in 1:max_doy
                line = readline(io)
                parts = split(line)
                if length(parts) >= 3
                    wert = parse(Float32, parts[2])
                    xlyman_obs[i] = wert
                else
                    error("Unexpected line format: $line")
                end
            end
        end
        println("      Lyman‑α read successfully.")
        println("      xlyman_obs first 5: ", xlyman_obs[1:5])
        println("      xlyman_obs last 5:  ", xlyman_obs[end-4:end])
    catch e
        println("      Warning: could not read Lyman‑α file, using default 4.5e11")
        fill!(xlyman_obs, 4.5f11)
    end

    # ------------------------------------------------------------------
    # 5. UT1 file – formatted, get unix_time for 10 May of this year
    # ------------------------------------------------------------------
    println("    Reading UT1 table: $(cfg.ut1_file)")
    unix_time = 0.0f0
    try
        open(cfg.ut1_file, "r") do io
            readline(io)   # skip header
            for i in 1:240   # Fortran reads 240 lines
                line = readline(io)
                parts = split(line)
                if length(parts) >= 2
                    b1 = parse(Float32, parts[1])
                    b2 = parse(Float32, parts[2])
                    if Float32(year) == b1
                        unix_time = b2
                    end
                end
            end
        end
        println("      unix_time = $unix_time")
    catch e
        println("      Warning: could not read UT1 file, using dummy 0.0")
        unix_time = 0.0f0
    end

    # ------------------------------------------------------------------
    # 6. CH4 file – formatted, compute scale_ch4
    # ------------------------------------------------------------------
    println("    Reading CH4 table: $(cfg.ch4_file)")
    ch4_akt = 0.0f0
    ch4_2008 = 0.0f0
    try
        open(cfg.ch4_file, "r") do io
            readline(io)   # skip header
            for i in 1:275   # Fortran reads 275 lines (1861‑2135)
                line = readline(io)
                parts = split(line)
                if length(parts) >= 2
                    a1 = parse(Float32, parts[1])
                    a2 = parse(Float32, parts[2])
                    if Float32(year - 5) == a1 
                        ch4_akt = a2
                    end
                    if a1 == 2008.0f0
                        ch4_2008 = a2
                    end
                end
            end
        end
    catch e
        println("      Warning: could not read CH4 file, using dummy values")
        ch4_akt = 0.0f0
        ch4_2008 = 0.0f0
    end

    fr_alpha = 0.95f0
    ch4_akt   = ch4_akt / 1000.0f0
    ch4_2008  = ch4_2008 / 1000.0f0
    h2oanteil_2008 = fr_alpha * ch4_2008 * 2.0f0
    h2oanteil_akt  = fr_alpha * ch4_akt * 2.0f0
    basis_h2o = 6.0f0 - h2oanteil_2008
    h2o_true  = basis_h2o + h2oanteil_akt
    scale_ch4 = 100.0f0 * h2o_true / 6.0f0   # in percent
    scale_ch4 = scale_ch4 * 0.01f0            # convert to fraction

    println("      scale_ch4 = $scale_ch4")

    println("  init_global! done.")
    return xlyman_obs, scale_ch4, unix_time
end

# ----------------------------------------------------------------------
# Read initial H2O pressure profile (sub_leseh2o_zprinit)
# ----------------------------------------------------------------------
"""
    read_h2o_init!(S::MIMASState, cfg::Config)

Fortran equivalent: `sub_leseh2o_zprinit` (called after global init).
Reads the file `H2OINIT` (from `cfg.h2oinit`) containing two arrays:
- zpress_init: pressure levels (km) of length kgitzpr (166 levels from 72 to 105 km)
- zpress_init_h2o: corresponding H2O mixing ratio (ppmv)
These are used later to interpolate the initial 3D H2O field onto the model grid.
If the file cannot be read, dummy values are used (linear pressure, constant 5 ppmv).
"""
function read_h2o_init!(S::MIMASState, cfg::Config)
    println("    Reading H2O initialisation file: $(cfg.h2oinit)")
    try
        open(cfg.h2oinit, "r") do io
            # The file contains two separate Fortran unformatted records,
            # each an array of Float32 of length kgitzpr.
            S.h2o.zpress_init    .= read_fortran_record(io, Float32, cfg.kgitzpr)
            S.h2o.zpress_init_h2o .= read_fortran_record(io, Float32, cfg.kgitzpr)
        end
        # Print a few values to verify the read was successful
        println("      zpress_init first 5: ", S.h2o.zpress_init[1:5])
        println("      zpress_init last 5:  ", S.h2o.zpress_init[end-4:end])
        println("      zpress_init_h2o first 5: ", S.h2o.zpress_init_h2o[1:5])
        println("      zpress_init_h2o last 5:  ", S.h2o.zpress_init_h2o[end-4:end])
    catch e
        println("      Warning: could not read H2O init file, using dummy values")
        # Dummy: linear pressure from 72 to 105 km, step 0.2 km, constant 5 ppmv
        for k in 1:cfg.kgitzpr
            S.h2o.zpress_init[k] = 72.0f0 + (k-1) * 0.2f0
        end
        fill!(S.h2o.zpress_init_h2o, 5.0f0)
        println("      dummy zpress_init first 5: ", S.h2o.zpress_init[1:5])
        println("      dummy zpress_init_h2o first 5: ", S.h2o.zpress_init_h2o[1:5])
    end
    println("    read_h2o_init! done.")
end


# ----------------------------------------------------------------------
# Read dynamics snapshot (sub_lesedyn) – full implementation with mean file corrections
# ----------------------------------------------------------------------
"""
    read_dynamics!(S::MIMASState, cfg::Config, ihemi::Int, ijahr_fixdyn::Int,
                   year::Int, month::Int, day::Int, hour::Int)

Fortran equivalent: `sub_lesedyn` (called at the start and then every hour).
Reads a dynamics snapshot (winds, temperature, density, pressure altitude) for the
given date/time and fills the `*1` arrays (`um1`, `tm1`, `dm1`, `vm1`, `wm1`, `zm1`)
in `S.uvwtd48`. The reading includes:
- Raw fields from the main file (dichte_lima, u_lima, v_lima, w_lima, t_lima)
- Pressure altitude from a separate file (zm)
- Corrections using mean files (1‑run and 6‑run climatologies)
- Temperature floor (102 K) and southern hemisphere mirroring
- Interpolation to cell centres to obtain all fields on the model grid
- Final smoothing of `wm1` (5 iterations of a 3×3 weighted average near the pole)
After this, the state contains the first of the two time levels needed for the time loop.
"""
function read_dynamics!(S, cfg, ihemi, ijahr_fixdyn, year, month, day, hour)
    println("    Reading dynamics for $(ijahr_fixdyn)-$(month)-$(day) $(hour):00")

    # ------------------------------------------------------------------
    # 1. Construct file name suffixes based on hemisphere
    # ------------------------------------------------------------------
    if ihemi == 1
        suffix_hem = ".hem_s_icegrid"
        suffix_zpr = ".zpr_s_icegrid"
    else
        suffix_hem = ".hem_n_icegrid"
        suffix_zpr = ".zpr_n_icegrid"
    end

    file_base = "$(lpad(ijahr_fixdyn,4,'0'))/$(lpad(month,2,'0'))/$(lpad(ijahr_fixdyn,4,'0'))-$(lpad(month,2,'0'))-$(lpad(day,2,'0')).$(lpad(hour,2,'0'))"
    file_zpr = joinpath(cfg.lima_pfad, file_base * suffix_zpr)
    file_hem = joinpath(cfg.lima_pfad, file_base * suffix_hem)

    # ------------------------------------------------------------------
    # 2. Read pressure altitude file (zm) – one record
    # ------------------------------------------------------------------
    println("      Reading pressure altitude: $file_zpr")
    try
        open(file_zpr, "r") do io
            S.uvwtd48.zm1 .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
        end
        # Print a sample to verify read
        println("        zm1[1,1,1] = ", S.uvwtd48.zm1[1,1,1])
        println("        zm1[53,120,163] = ", S.uvwtd48.zm1[end,end,end])
    catch e
        error("Failed to read pressure altitude file: $file_zpr\n", e)
    end

    # ------------------------------------------------------------------
    # 3. Read dynamics file (hem) – 5 separate records
    # ------------------------------------------------------------------
    println("      Reading dynamics: $file_hem")
    dichte_lima = Array{Float32}(undef, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
    u_lima      = similar(dichte_lima)
    v_lima      = similar(dichte_lima)
    w_lima      = similar(dichte_lima)
    t_lima      = similar(dichte_lima)

    try
        open(file_hem, "r") do io
            dichte_lima .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
            u_lima      .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
            v_lima      .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
            w_lima      .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
            t_lima      .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
        end
        # Print a sample from one field to verify
        println("        u_lima[1,1,1] = ", u_lima[1,1,1])
        println("        t_lima[1,1,1] = ", t_lima[1,1,1])
    catch e
        error("Failed to read dynamics file: $file_hem\n", e)
    end

    # ------------------------------------------------------------------
    # 4. Read mean files (5 records each) – these are used for corrections
    # ------------------------------------------------------------------
    # Helper to read a mean file containing 5 separate records of shape (31, kgitneu, nbneu)
    function read_mean_file(yr, mo, suffix)
        file = joinpath(cfg.lima_pfad, "$(lpad(yr,4,'0'))/$(lpad(mo,2,'0'))/uvwtd_mean_$(lpad(yr,4,'0'))_$(lpad(mo,2,'0'))$(suffix)")
        arrs = Vector{Array{Float32,3}}(undef, 5)
        try
            open(file, "r") do io
                arrs[1] = read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
                arrs[2] = read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
                arrs[3] = read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
                arrs[4] = read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
                arrs[5] = read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
            end
        catch e
            error("Failed to read mean file: $file\n", e)
        end
        return (arrs[1], arrs[2], arrs[3], arrs[4], arrs[5])
    end

    # Read 1‑run fix for fixed year (ijahr_fixdyn)
    println("        Reading 1‑run mean files for $(ijahr_fixdyn)-$(month)")
    dm_out_1run_fix, um_out_1run_fix, vm_out_1run_fix, wm_out_1run_fix, tm_out_1run_fix =
        read_mean_file(ijahr_fixdyn, month, ".bin")

    # Read 6‑run fix for fixed year
    println("        Reading 6‑run mean files for $(ijahr_fixdyn)-$(month)")
    dm_out_6run_fix, um_out_6run_fix, vm_out_6run_fix, wm_out_6run_fix, tm_out_6run_fix =
        read_mean_file(ijahr_fixdyn, month, ".binzz")

    # Read 6‑run for current year (if different)
    if year == ijahr_fixdyn
        # Reuse the arrays to avoid reading again
        dm_out_6run = dm_out_6run_fix
        um_out_6run = um_out_6run_fix
        vm_out_6run = vm_out_6run_fix
        wm_out_6run = wm_out_6run_fix
        tm_out_6run = tm_out_6run_fix
        println("        Using same 6‑run means for current year (year == ijahr_fixdyn)")
    else
        println("        Reading 6‑run mean files for $(year)-$(month)")
        dm_out_6run, um_out_6run, vm_out_6run, wm_out_6run, tm_out_6run =
            read_mean_file(year, month, ".binzz")
    end

    # ------------------------------------------------------------------
    # 5. Read pressure altitude mean files (single record each)
    # ------------------------------------------------------------------
    function read_zpress_mean(yr, mo, suffix)
        file = joinpath(cfg.lima_pfad, "$(lpad(yr,4,'0'))/$(lpad(mo,2,'0'))/zpress_mean_$(lpad(yr,4,'0'))_$(lpad(mo,2,'0'))$(suffix)")
        try
            open(file, "r") do io
                return read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
            end
        catch e
            error("Failed to read pressure altitude mean file: $file\n", e)
        end
    end

    println("        Reading pressure altitude means for $(ijahr_fixdyn)-$(month) (1‑run)")
    zm_out_1run_fix = read_zpress_mean(ijahr_fixdyn, month, ".bin_zpr")
    println("        Reading pressure altitude means for $(ijahr_fixdyn)-$(month) (6‑run)")
    zm_out_6run_fix = read_zpress_mean(ijahr_fixdyn, month, ".bin_zprzz")
    println("        Reading pressure altitude means for $(year)-$(month) (6‑run)")
    zm_out_6run     = read_zpress_mean(year, month, ".bin_zprzz")

    # ------------------------------------------------------------------
    # 6. Correct the fields using the mean files
    # ------------------------------------------------------------------
    itag = day   # day of month (1‑based)

    println("        Applying corrections with day index itag = ", itag)
    for k in 1:cfg.kgitneu
        for i in 1:cfg.igitneu
            for j in 1:cfg.nbneu
                # Winds: use 6run_fix (permanently fixed year)
                u_lima[j,i,k] = u_lima[j,i,k] - um_out_1run_fix[itag,k,j] + um_out_6run_fix[itag,k,j]
                v_lima[j,i,k] = v_lima[j,i,k] - vm_out_1run_fix[itag,k,j] + vm_out_6run_fix[itag,k,j]
                w_lima[j,i,k] = w_lima[j,i,k] - wm_out_1run_fix[itag,k,j] + wm_out_6run_fix[itag,k,j]

                # Density, temperature, pressure altitude: use 6run (current year)
                dichte_lima[j,i,k] = dichte_lima[j,i,k] - dm_out_1run_fix[itag,k,j] + dm_out_6run[itag,k,j]
                t_lima[j,i,k]      = t_lima[j,i,k]      - tm_out_1run_fix[itag,k,j] + tm_out_6run[itag,k,j]
                S.uvwtd48.zm1[j,i,k] = S.uvwtd48.zm1[j,i,k] - zm_out_1run_fix[itag,k,j] + zm_out_6run[itag,k,j]
            end
        end
    end

    # ------------------------------------------------------------------
    # 7. Temperature floor (Fortran: if t < 102 K, set to 102 K)
    # ------------------------------------------------------------------
    for k in 1:cfg.kgitneu, i in 1:cfg.igitneu, j in 1:cfg.nbneu
        if t_lima[j,i,k] < 102.0f0
            t_lima[j,i,k] = 102.0f0
        end
    end

    # ------------------------------------------------------------------
    # 8. Mirror winds for southern hemisphere (ihemi == 1)
    # ------------------------------------------------------------------
    if ihemi == 1
        v_lima .= -v_lima
        println("        Southern hemisphere: v_lima mirrored (negated).")
    end

    # ------------------------------------------------------------------
    # 9. Interpolate to cell centers to get um, tm, dm, vm, wm
    # ------------------------------------------------------------------
    # Create temporary arrays on the shifted grid
    w_shift = Array{Float32}(undef, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
    tm_shift = similar(w_shift)
    v_shift = similar(w_shift)
    dichte_shift = similar(w_shift)

    for k in 1:cfg.kgitneu
        for j in 1:cfg.nbneu
            for i in 1:cfg.igitneu-1
                w_shift[j,i,k]   = 0.5f0 * (w_lima[j,i,k] + w_lima[j,i+1,k])
                tm_shift[j,i,k]  = 0.5f0 * (t_lima[j,i,k] + t_lima[j,i+1,k])
                v_shift[j,i,k]   = 0.5f0 * (v_lima[j,i,k] + v_lima[j,i+1,k])
                dichte_shift[j,i,k] = 0.5f0 * (dichte_lima[j,i,k] + dichte_lima[j,i+1,k])
            end
            # i = igitneu (wrap around)
            i = cfg.igitneu
            w_shift[j,i,k]   = 0.5f0 * (w_lima[j,i,k] + w_lima[j,1,k])
            tm_shift[j,i,k]  = 0.5f0 * (t_lima[j,i,k] + t_lima[j,1,k])
            v_shift[j,i,k]   = 0.5f0 * (v_lima[j,i,k] + v_lima[j,1,k])
            dichte_shift[j,i,k] = 0.5f0 * (dichte_lima[j,i,k] + dichte_lima[j,1,k])
        end
    end

    # Compute vm (on nbneu1 grid)
    vm1 = Array{Float32}(undef, cfg.nbneu1, cfg.igitneu, cfg.kgitneu)
    for k in 1:cfg.kgitneu
        for i in 1:cfg.igitneu
            for j in 2:cfg.nbneu
                vm1[j,i,k] = 0.5f0 * (v_shift[j-1,i,k] + v_shift[j,i,k])
            end
            vm1[1,i,k] = vm1[2,i,k]          # j=1 set to j=2
            vm1[cfg.nbneu1,i,k] = 0.0f0      # top boundary zero
        end
    end

    # Compute wm (on kgitneu1 grid)
    wm1 = Array{Float32}(undef, cfg.nbneu, cfg.igitneu, cfg.kgitneu1)
    for i in 1:cfg.igitneu
        for j in 1:cfg.nbneu
            for k in 2:cfg.kgitneu
                wm1[j,i,k] = 0.5f0 * (w_shift[j,i,k-1] + w_shift[j,i,k])
            end
            wm1[j,i,1] = wm1[j,i,2]                     # bottom duplicate
            wm1[j,i,cfg.kgitneu1] = wm1[j,i,cfg.kgitneu] # top duplicate
        end
    end

    # Compute pressure dm from density and temperature
    rstern = 8.31432f3   # 8.31432 * 1000 (universal gas constant in J/(kmol·K) * 1000? Fortran uses this.)
    dm1 = Array{Float32}(undef, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
    for k in 1:cfg.kgitneu
        for i in 1:cfg.igitneu
            for j in 1:cfg.nbneu
                dm1[j,i,k] = rstern * dichte_shift[j,i,k] * tm_shift[j,i,k] / 28.96f0
            end
        end
    end

    # ------------------------------------------------------------------
    # 10. Store results into S.uvwtd48
    # ------------------------------------------------------------------
    S.uvwtd48.um1 .= u_lima   # u_lima is on original grid (no shift)
    S.uvwtd48.tm1 .= tm_shift
    S.uvwtd48.dm1 .= dm1
    S.uvwtd48.vm1 .= vm1
    S.uvwtd48.wm1 .= wm1
    # zm1 already updated in correction step

    # Print a sample to verify final fields
    println("        Final um1[1,1,1] = ", S.uvwtd48.um1[1,1,1])
    println("        Final tm1[1,1,1] = ", S.uvwtd48.tm1[1,1,1])
    println("        Final dm1[1,1,1] = ", S.uvwtd48.dm1[1,1,1])
    println("        Final vm1[1,1,1] = ", S.uvwtd48.vm1[1,1,1])
    println("        Final wm1[1,1,1] = ", S.uvwtd48.wm1[1,1,1])

    # ------------------------------------------------------------------
    # 11. Apply smoothing to wm1 (5 iterations + top adjustment)
    # ------------------------------------------------------------------
    println("        Applying 5‑step smoothing to wm1 near the pole...")
    smooth_wm!(S.uvwtd48.wm1, cfg)
    println("        wm1[1,1,1] after smoothing = ", S.uvwtd48.wm1[1,1,1])

    println("    read_dynamics! done.")
    return nothing
end


# ----------------------------------------------------------------------
# Smoothing of vertical wind field (helper for read_dynamics!)
# ----------------------------------------------------------------------
"""
    smooth_wm!(wm::Array{Float32,3}, cfg::Config)

Apply the 5‑iteration horizontal smoothing to `wm` near the pole,
followed by the linear combination for the highest latitudes.
This replicates the Fortran code that smooths `wm1` after reading.
The smoothing uses a weighted 3×3 stencil (centre weight 27, neighbours 1, total 35)
over the region j = nbneu-9 .. nbneu-4 (i.e., latitudes near the North Pole).
After 5 iterations, a linear combination is applied to the topmost latitudes
(j = 50..53) using values from j = 48 and 49.
"""
function smooth_wm!(wm, cfg)
    nbneu = cfg.nbneu
    igitneu = cfg.igitneu
    kgitneu1 = cfg.kgitneu1
    workw = similar(wm)

    for iter in 1:5
        workw .= wm
        for k in 2:kgitneu1-1
            for i in 1:igitneu
                i0 = i == 1 ? igitneu : i-1
                i2 = i == igitneu ? 1 : i+1
                for j in nbneu-9:nbneu-4
                    wm[j,i,k] = (
                        workw[j-1, i0, k] + workw[j, i0, k] + workw[j+1, i0, k] +
                        workw[j-1, i,  k] + 27.0f0 * workw[j, i, k] + workw[j+1, i, k] +
                        workw[j-1, i2, k] + workw[j, i2, k] + workw[j+1, i2, k]
                    ) / 35.0f0
                end
            end
        end
    end

    # Top‑level adjustment (j=50..53)
    for k in 1:kgitneu1
        for i in 1:igitneu
            a1 = 0.4f0 * wm[49, i, k]
            a2 = 0.4f0 * wm[48, i, k]
            wm[50, i, k] = 0.2f0 * wm[50, i, k] + a1 + a2
            wm[51, i, k] = 0.2f0 * wm[51, i, k] + a1 + a2
            wm[52, i, k] = 0.2f0 * wm[52, i, k] + a1 + a2
            wm[53, i, k] = 0.2f0 * wm[53, i, k] + a1 + a2
        end
    end
    nothing
end

# ----------------------------------------------------------------------
# Writing restart file (Fortran unformatted output)
# ----------------------------------------------------------------------
# These functions write data in the same Fortran unformatted format that
# read_fortran_record expects, using big‑endian byte order and 4‑byte record
# headers/trailers. They are used to write the restart file after a cold start.

"""
    write_fortran_record(io::IO, data::Array{T,N}) where {T,N}

Write an array `data` to `io` as a Fortran unformatted record.
- Writes a 4‑byte header containing the total number of bytes in the record (big‑endian).
- Writes the array data, converting to big‑endian if the machine is little‑endian.
- Writes a 4‑byte trailer (same as header).
This matches the format produced by Fortran with `-convert big_endian`.
"""
function write_fortran_record(io::IO, data::Array{T,N}) where {T,N}
    # Compute total bytes to be written
    nbytes = sizeof(T) * length(data)

    # Write header (big‑endian)
    write(io, ntoh(Int32(nbytes)))

    # Convert data to big‑endian if necessary
    if ENDIAN_BOM == 0x04030201   # little‑endian machine (most common)
        # Swap bytes for each element (bswap works for Float32 as well)
        data_be = bswap.(data)
        write(io, data_be)
    else
        # Machine is big‑endian, write as is
        write(io, data)
    end

    # Write trailer (same as header)
    write(io, ntoh(Int32(nbytes)))
end


"""
    write_restart!(S::MIMASState, cfg::Config)

Fortran equivalent: the block of code after `sub_init_staub` that writes the
restart file `CRINIT` when `nstart == 0`.
Writes a restart file containing:
- An integer `i_std = 0` (written as a record containing a single Int32)
- The four particle arrays `xfeld`, `yfeld`, `zfeld`, `rfeld` (each as a separate record)
The file is written in Fortran unformatted format (big‑endian with record markers).
If the write fails, it retries every 2 seconds (as in the original Fortran code).
"""
function write_restart!(S, cfg)
    if cfg.nstart != 0
        return   # only for cold start
    end

    i_std = 0
    println("    Writing restart file: $(cfg.crinit)")
    local io
    success = false
    while !success
        try
            io = open(cfg.crinit, "w")
            # Write integer i_std as a record (Fortran writes a single integer)
            write_fortran_record(io, [Int32(i_std)])   # wrapped in array for uniform handling

            # Write particle arrays
            write_fortran_record(io, S.eis.xfeld)
            write_fortran_record(io, S.eis.yfeld)
            write_fortran_record(io, S.eis.zfeld)
            write_fortran_record(io, S.eis.rfeld)

            close(io)
            success = true
        catch e
            println("      Error writing restart file: $e. Retrying in 2 seconds...")
            sleep(2)
            # Continue loop
        end
    end
    # Print confirmation and file size (if file exists)
    if isfile(cfg.crinit)
        filesize_mb = round(filesize(cfg.crinit) / 1e6, digits=2)
        println("    Restart file written successfully (size: $(filesize_mb) MB).")
    else
        println("    Restart file written successfully.")
    end
end

# ----------------------------------------------------------------------
# Read dynamics snapshot into arbitrary target arrays (for *2)
# ----------------------------------------------------------------------
"""
    read_dynamics_into!(um, tm, dm, vm, wm, zm, S, cfg, ihemi, ijahr_fixdyn,
                        year, month, day, hour)

Read a dynamics snapshot for the given date/time and fill the provided target arrays
(um, tm, dm, vm, wm, zm) with the processed fields. The arrays must have the correct
dimensions. This is used to read the next hour's data into the *2 arrays.
"""
function read_dynamics_into!(um, tm, dm, vm, wm, zm, S, cfg, ihemi, ijahr_fixdyn,
                             year, month, day, hour)
    println("      Reading dynamics into target arrays for $(ijahr_fixdyn)-$(month)-$(day) $(hour):00")

    # Construct file names (same as in read_dynamics!)
    if ihemi == 1
        suffix_hem = ".hem_s_icegrid"
        suffix_zpr = ".zpr_s_icegrid"
    else
        suffix_hem = ".hem_n_icegrid"
        suffix_zpr = ".zpr_n_icegrid"
    end

    file_base = "$(lpad(ijahr_fixdyn,4,'0'))/$(lpad(month,2,'0'))/$(lpad(ijahr_fixdyn,4,'0'))-$(lpad(month,2,'0'))-$(lpad(day,2,'0')).$(lpad(hour,2,'0'))"
    file_zpr = joinpath(cfg.lima_pfad, file_base * suffix_zpr)
    file_hem = joinpath(cfg.lima_pfad, file_base * suffix_hem)

    # ------------------------------------------------------------------
    # 1. Read pressure altitude file (zm) – one record
    # ------------------------------------------------------------------
    println("        Reading pressure altitude: $file_zpr")
    try
        open(file_zpr, "r") do io
            zm .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
        end
    catch e
        error("Failed to read pressure altitude file: $file_zpr\n", e)
    end

    # ------------------------------------------------------------------
    # 2. Read dynamics file (hem) – 5 separate records
    # ------------------------------------------------------------------
    println("        Reading dynamics: $file_hem")
    dichte_lima = Array{Float32}(undef, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
    u_lima      = similar(dichte_lima)
    v_lima      = similar(dichte_lima)
    w_lima      = similar(dichte_lima)
    t_lima      = similar(dichte_lima)

    try
        open(file_hem, "r") do io
            dichte_lima .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
            u_lima      .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
            v_lima      .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
            w_lima      .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
            t_lima      .= read_fortran_record(io, Float32, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
        end
    catch e
        error("Failed to read dynamics file: $file_hem\n", e)
    end

    # ------------------------------------------------------------------
    # 3. Read mean files (same as in read_dynamics!)
    # ------------------------------------------------------------------
    # Helper to read a mean file containing 5 separate records of shape (31, kgitneu, nbneu)
    function read_mean_file(yr, mo, suffix)
        file = joinpath(cfg.lima_pfad, "$(lpad(yr,4,'0'))/$(lpad(mo,2,'0'))/uvwtd_mean_$(lpad(yr,4,'0'))_$(lpad(mo,2,'0'))$(suffix)")
        arrs = Vector{Array{Float32,3}}(undef, 5)
        try
            open(file, "r") do io
                arrs[1] = read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
                arrs[2] = read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
                arrs[3] = read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
                arrs[4] = read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
                arrs[5] = read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
            end
        catch e
            error("Failed to read mean file: $file\n", e)
        end
        return (arrs[1], arrs[2], arrs[3], arrs[4], arrs[5])
    end

    # Read 1‑run fix for fixed year (ijahr_fixdyn)
    println("        Reading 1‑run mean files for $(ijahr_fixdyn)-$(month)")
    dm_out_1run_fix, um_out_1run_fix, vm_out_1run_fix, wm_out_1run_fix, tm_out_1run_fix =
        read_mean_file(ijahr_fixdyn, month, ".bin")

    # Read 6‑run fix for fixed year
    println("        Reading 6‑run mean files for $(ijahr_fixdyn)-$(month)")
    dm_out_6run_fix, um_out_6run_fix, vm_out_6run_fix, wm_out_6run_fix, tm_out_6run_fix =
        read_mean_file(ijahr_fixdyn, month, ".binzz")

    # Read 6‑run for current year (if different)
    if year == ijahr_fixdyn
        dm_out_6run = dm_out_6run_fix
        um_out_6run = um_out_6run_fix
        vm_out_6run = vm_out_6run_fix
        wm_out_6run = wm_out_6run_fix
        tm_out_6run = tm_out_6run_fix
        println("        Using same 6‑run means for current year (year == ijahr_fixdyn)")
    else
        println("        Reading 6‑run mean files for $(year)-$(month)")
        dm_out_6run, um_out_6run, vm_out_6run, wm_out_6run, tm_out_6run =
            read_mean_file(year, month, ".binzz")
    end

    # ------------------------------------------------------------------
    # 4. Read pressure altitude mean files (single record each)
    # ------------------------------------------------------------------
    function read_zpress_mean(yr, mo, suffix)
        file = joinpath(cfg.lima_pfad, "$(lpad(yr,4,'0'))/$(lpad(mo,2,'0'))/zpress_mean_$(lpad(yr,4,'0'))_$(lpad(mo,2,'0'))$(suffix)")
        try
            open(file, "r") do io
                return read_fortran_record(io, Float32, 31, cfg.kgitneu, cfg.nbneu)
            end
        catch e
            error("Failed to read pressure altitude mean file: $file\n", e)
        end
    end

    println("        Reading pressure altitude means for $(ijahr_fixdyn)-$(month) (1‑run)")
    zm_out_1run_fix = read_zpress_mean(ijahr_fixdyn, month, ".bin_zpr")
    println("        Reading pressure altitude means for $(ijahr_fixdyn)-$(month) (6‑run)")
    zm_out_6run_fix = read_zpress_mean(ijahr_fixdyn, month, ".bin_zprzz")
    println("        Reading pressure altitude means for $(year)-$(month) (6‑run)")
    zm_out_6run     = read_zpress_mean(year, month, ".bin_zprzz")

    # ------------------------------------------------------------------
    # 5. Correct the fields using the mean files
    # ------------------------------------------------------------------
    itag = day   # day of month (1‑based)

    println("        Applying corrections with day index itag = ", itag)
    for k in 1:cfg.kgitneu
        for i in 1:cfg.igitneu
            for j in 1:cfg.nbneu
                # Winds: use 6run_fix
                u_lima[j,i,k] = u_lima[j,i,k] - um_out_1run_fix[itag,k,j] + um_out_6run_fix[itag,k,j]
                v_lima[j,i,k] = v_lima[j,i,k] - vm_out_1run_fix[itag,k,j] + vm_out_6run_fix[itag,k,j]
                w_lima[j,i,k] = w_lima[j,i,k] - wm_out_1run_fix[itag,k,j] + wm_out_6run_fix[itag,k,j]

                # Density, temperature, pressure altitude: use 6run (current year)
                dichte_lima[j,i,k] = dichte_lima[j,i,k] - dm_out_1run_fix[itag,k,j] + dm_out_6run[itag,k,j]
                t_lima[j,i,k]      = t_lima[j,i,k]      - tm_out_1run_fix[itag,k,j] + tm_out_6run[itag,k,j]
                # Note: we do NOT modify S.uvwtd48.zm1 here; instead we use the passed zm array,
                # which has already been read. The original Fortran code updates zm1 in place.
                # Here we update the local zm array (the target).
                zm[j,i,k] = zm[j,i,k] - zm_out_1run_fix[itag,k,j] + zm_out_6run[itag,k,j]
            end
        end
    end

    # ------------------------------------------------------------------
    # 6. Temperature floor
    # ------------------------------------------------------------------
    for k in 1:cfg.kgitneu, i in 1:cfg.igitneu, j in 1:cfg.nbneu
        if t_lima[j,i,k] < 102.0f0
            t_lima[j,i,k] = 102.0f0
        end
    end

    # ------------------------------------------------------------------
    # 7. Mirror winds for southern hemisphere
    # ------------------------------------------------------------------
    if ihemi == 1
        v_lima .= -v_lima
        println("        Southern hemisphere: v_lima mirrored (negated).")
    end

    # ------------------------------------------------------------------
    # 8. Interpolate to cell centers to get final fields
    # ------------------------------------------------------------------
    # Create temporary arrays on the shifted grid
    w_shift = Array{Float32}(undef, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
    tm_shift = similar(w_shift)
    v_shift = similar(w_shift)
    dichte_shift = similar(w_shift)

    for k in 1:cfg.kgitneu
        for j in 1:cfg.nbneu
            for i in 1:cfg.igitneu-1
                w_shift[j,i,k]   = 0.5f0 * (w_lima[j,i,k] + w_lima[j,i+1,k])
                tm_shift[j,i,k]  = 0.5f0 * (t_lima[j,i,k] + t_lima[j,i+1,k])
                v_shift[j,i,k]   = 0.5f0 * (v_lima[j,i,k] + v_lima[j,i+1,k])
                dichte_shift[j,i,k] = 0.5f0 * (dichte_lima[j,i,k] + dichte_lima[j,i+1,k])
            end
            i = cfg.igitneu
            w_shift[j,i,k]   = 0.5f0 * (w_lima[j,i,k] + w_lima[j,1,k])
            tm_shift[j,i,k]  = 0.5f0 * (t_lima[j,i,k] + t_lima[j,1,k])
            v_shift[j,i,k]   = 0.5f0 * (v_lima[j,i,k] + v_lima[j,1,k])
            dichte_shift[j,i,k] = 0.5f0 * (dichte_lima[j,i,k] + dichte_lima[j,1,k])
        end
    end

    # Compute vm (on nbneu1 grid)
    vm1 = Array{Float32}(undef, cfg.nbneu1, cfg.igitneu, cfg.kgitneu)
    for k in 1:cfg.kgitneu
        for i in 1:cfg.igitneu
            for j in 2:cfg.nbneu
                vm1[j,i,k] = 0.5f0 * (v_shift[j-1,i,k] + v_shift[j,i,k])
            end
            vm1[1,i,k] = vm1[2,i,k]          # j=1 set to j=2
            vm1[cfg.nbneu1,i,k] = 0.0f0      # top boundary zero
        end
    end

    # Compute wm (on kgitneu1 grid)
    wm1 = Array{Float32}(undef, cfg.nbneu, cfg.igitneu, cfg.kgitneu1)
    for i in 1:cfg.igitneu
        for j in 1:cfg.nbneu
            for k in 2:cfg.kgitneu
                wm1[j,i,k] = 0.5f0 * (w_shift[j,i,k-1] + w_shift[j,i,k])
            end
            wm1[j,i,1] = wm1[j,i,2]                     # bottom duplicate
            wm1[j,i,cfg.kgitneu1] = wm1[j,i,cfg.kgitneu] # top duplicate
        end
    end

    # Compute pressure dm from density and temperature
    rstern = 8.31432f3
    dm1 = Array{Float32}(undef, cfg.nbneu, cfg.igitneu, cfg.kgitneu)
    for k in 1:cfg.kgitneu
        for i in 1:cfg.igitneu
            for j in 1:cfg.nbneu
                dm1[j,i,k] = rstern * dichte_shift[j,i,k] * tm_shift[j,i,k] / 28.96f0
            end
        end
    end

    # ------------------------------------------------------------------
    # 9. Store results into target arrays
    # ------------------------------------------------------------------
    um .= u_lima          # u_lima is on original grid (no shift)
    tm .= tm_shift
    dm .= dm1
    vm .= vm1
    wm .= wm1
    # zm already updated in correction step

    println("        read_dynamics_into! done.")
end

# ----------------------------------------------------------------------
# Output: schreibe_beta0! – fill netCDF buffers
# ----------------------------------------------------------------------
"""
    schreibe_beta0!(S, cfg, unix_time)

Fortran equivalent: `schreibe_beta0`.
Accumulates particle data and dynamics fields into the netCDF buffers
(S.netcdf.aq, bq, ...) for the current output time. Increments S.tu.
"""
function schreibe_beta0!(S, cfg, unix_time)
    S.tu += 1
    S.times[S.tu] = unix_time

    nbneu = cfg.nbneu
    igitneu = cfg.igitneu
    kgitneu = cfg.kgitneu
    ntrac = cfg.ntrac

    # Local temporary arrays for accumulation
    xneis = zeros(Float32, nbneu, igitneu, kgitneu)
    xnstaub = zeros(Float32, nbneu, igitneu, kgitneu)
    eisradius = zeros(Float32, nbneu, igitneu, kgitneu)
    eisradiuseff = zeros(Float32, nbneu, igitneu, kgitneu)
    reff_zaehler = zeros(Float32, nbneu, igitneu, kgitneu)
    reff_nenner = zeros(Float32, nbneu, igitneu, kgitneu)
    betaorg = zeros(Float32, nbneu, igitneu, kgitneu)
    ext_coef126 = zeros(Float32, nbneu, igitneu, kgitneu)
    ext_coef200 = zeros(Float32, nbneu, igitneu, kgitneu)
    ext_coef532 = zeros(Float32, nbneu, igitneu, kgitneu)
    ext_coef1000 = zeros(Float32, nbneu, igitneu, kgitneu)
    ext_coef3000 = zeros(Float32, nbneu, igitneu, kgitneu)
    densice = zeros(Float32, nbneu, igitneu, kgitneu)

    # Counters for diagnostics (optional)
    neis = 0
    nstaub = 0

    xpi = pi
    rhoicecgs = rhoice * Float32(1e-3)   # convert kg/m³ to g/cm³

    # First pass: accumulate particle contributions
    for n in 1:ntrac
        i = floor(Int, S.eis.xfeld[n])
        j = floor(Int, S.eis.yfeld[n])
        k = floor(Int, S.eis.zfeld[n])

        i = clamp(i, 1, igitneu)
        j = clamp(j, 1, nbneu)
        k = clamp(k, 1, kgitneu)

        if S.eis.rfeld[n] > S.eis.rinit[n]
            # ice particle
            neis += 1
            xneis[j,i,k] += 1.0f0
            eisradius[j,i,k] += S.eis.rfeld[n]

            radi2 = S.eis.rfeld[n]^2
            radi3 = S.eis.rfeld[n] * radi2
            reff_zaehler[j,i,k] += radi3
            reff_nenner[j,i,k]  += radi2

            volume = (4.0f0/3.0f0) * xpi * radi3 * Float32(1e-21)   # cm³
            densice[j,i,k] += volume * rhoicecgs * Float32(1e15)      # g/km³

            icross = floor(Int, S.eis.rfeld[n] * 10.0f0 + 1.1f0)
            icross = clamp(icross, 1, 1500)

            betaorg[j,i,k] += Float32(1e16) * S.consts.crosssec532[icross]
            ext_coef126[j,i,k] += Float32(1e6) * S.consts.ext126[icross]
            ext_coef200[j,i,k] += Float32(1e6) * S.consts.ext200[icross]
            ext_coef532[j,i,k] += Float32(1e6) * S.consts.ext532[icross]
            ext_coef1000[j,i,k] += Float32(1e6) * S.consts.ext1000[icross]
            ext_coef3000[j,i,k] += Float32(1e6) * S.consts.ext3000[icross]
        else
            # dust particle
            nstaub += 1
            xnstaub[j,i,k] += 1.0f0
        end
    end

    # Compute mean radius where there are ice particles
    for k in 1:kgitneu
        for i in 1:igitneu
            for j in 1:nbneu
                if xneis[j,i,k] >= 1.0f0
                    eisradius[j,i,k] /= xneis[j,i,k]
                end
            end
        end
    end

    # Compute effective radius
    for k in 1:kgitneu
        for i in 1:igitneu
            for j in 1:nbneu
                if reff_nenner[j,i,k] >= 1.0f0
                    eisradiuseff[j,i,k] = reff_zaehler[j,i,k] / reff_nenner[j,i,k]
                end
            end
        end
    end

    # Latitude weighting factor (bgewicht)
    bgewicht = zeros(Float32, nbneu)
    for j in 1:nbneu
        lat = 37.5f0 + (j-1)
        bgewicht[j] = 1.0f0 / cosd(lat)
    end

    # Apply scaling and convert to physical units
    volufak = 1.0f0   # as in Fortran (commented out)
    volufak *= igitneu

    for k in 1:kgitneu
        for i in 1:igitneu
            for j in 1:nbneu
                scale = bgewicht[j] * volufak / 100.0f0
                xneis[j,i,k] *= scale
                xnstaub[j,i,k] *= scale
                betaorg[j,i,k] *= scale
                ext_coef126[j,i,k] *= scale
                ext_coef200[j,i,k] *= scale
                ext_coef532[j,i,k] *= scale
                ext_coef1000[j,i,k] *= scale
                ext_coef3000[j,i,k] *= scale
                densice[j,i,k] *= scale
            end
        end
    end

    # Store into netCDF buffers (note the order: (k,i,j,tu) as in Fortran)
    for j in 1:nbneu
        for i in 1:igitneu
            for k in 1:kgitneu
                S.netcdf.aq[k,i,j,S.tu] = betaorg[j,i,k]
                S.netcdf.bq[k,i,j,S.tu] = densice[j,i,k]
                S.netcdf.cq[k,i,j,S.tu] = eisradius[j,i,k]
                S.netcdf.dq[k,i,j,S.tu] = xneis[j,i,k]
                S.netcdf.eq[k,i,j,S.tu] = xnstaub[j,i,k]
                S.netcdf.fq[k,i,j,S.tu] = eisradiuseff[j,i,k]
                S.netcdf.gq[k,i,j,S.tu] = S.h2o.hm[j,i,k]
                S.netcdf.hq[k,i,j,S.tu] = S.uvwtd.um[j,i,k]
                S.netcdf.iq[k,i,j,S.tu] = S.uvwtd.vm[j,i,k]
                S.netcdf.jq[k,i,j,S.tu] = S.uvwtd.tm[j,i,k]
                S.netcdf.kq[k,i,j,S.tu] = S.uvwtd48.zm1[j,i,k]
                S.netcdf.lq[k,i,j,S.tu] = ext_coef126[j,i,k]
                S.netcdf.mq[k,i,j,S.tu] = ext_coef200[j,i,k]
                S.netcdf.nq[k,i,j,S.tu] = ext_coef532[j,i,k]
                S.netcdf.oq[k,i,j,S.tu] = ext_coef1000[j,i,k]
                S.netcdf.pq[k,i,j,S.tu] = ext_coef3000[j,i,k]
            end
        end
    end

    println("    Output at hour ", S.tu, " (unix_time = ", unix_time, ")")
    nothing
end

# ----------------------------------------------------------------------
# Write netCDF output (writegrid)
# ----------------------------------------------------------------------
"""
    write_netcdf!(S, cfg)

Fortran equivalent: `writegrid`.
Writes the netCDF file with all accumulated data.
"""
# ----------------------------------------------------------------------
# Write netCDF output (writegrid)
# ----------------------------------------------------------------------

function write_netcdf!(S, cfg)
    outfile = joinpath(cfg.output_pfad, string(cfg.start_year) * ".nc")
    println("Writing netCDF file: ", outfile)
    if isfile(outfile)
        rm(outfile)
    end

    # Define coordinate variables and attributes (nccreate also defines dimensions).
    NetCDF.nccreate(outfile, "lat", "lat", cfg.nbneu;
                    t = NetCDF.NC_FLOAT,
                    atts = Dict("units" => "degree_north",
                                "standard_name" => "latitude",
                                "axis" => "X"))
    NetCDF.nccreate(outfile, "lon", "lon", cfg.igitneu;
                    t = NetCDF.NC_FLOAT,
                    atts = Dict("units" => "degree_east",
                                "standard_name" => "longitude",
                                "axis" => "Y"))
    NetCDF.nccreate(outfile, "height", "height", cfg.kgitneu;
                    t = NetCDF.NC_FLOAT,
                    atts = Dict("units" => "km",
                                "standard_name" => "geometric height",
                                "axis" => "Z"))
    NetCDF.nccreate(outfile, "time", "time", S.tu;
                    t = NetCDF.NC_FLOAT,
                    atts = Dict("units" => "seconds since 1970-01-01 00:00:00",
                                "axis" => "T"))

    # Define data variables
    var_names = ["beta532", "densice", "r_ice", "n_ice", "n_dust", "r_ice_eff",
                 "H2O", "zon_win", "mer_win", "tem", "press_alt",
                 "ext_126", "ext_200", "ext_532", "ext_1000", "ext_3000"]
    data_arrays = [S.netcdf.aq, S.netcdf.bq, S.netcdf.cq, S.netcdf.dq, S.netcdf.eq,
                   S.netcdf.fq, S.netcdf.gq, S.netcdf.hq, S.netcdf.iq, S.netcdf.jq,
                   S.netcdf.kq, S.netcdf.lq, S.netcdf.mq, S.netcdf.nq, S.netcdf.oq,
                   S.netcdf.pq]
    units = [
        "10E-10m^-1sr^-1",
        "g/km^3",
        "nm",
        "cm^-3",
        "cm^-3",
        "nm",
        "ppmv",
        "m/s",
        "m/s",
        "K",
        "km",
        "m^-1",
        "m^-1",
        "m^-1",
        "m^-1",
        "m^-1"
    ]

    for (idx, name) in enumerate(var_names)
        NetCDF.nccreate(outfile, name,
                        "height", cfg.kgitneu,
                        "lon", cfg.igitneu,
                        "lat", cfg.nbneu,
                        "time", S.tu;
                        t = NetCDF.NC_FLOAT,
                        atts = Dict("units" => units[idx]))
    end

    # Write coordinate data
    lats = [37.5f0 + (j-1) for j in 1:cfg.nbneu]
    lons = [Float32((i-1)*3) for i in 1:cfg.igitneu]
    NetCDF.ncwrite(lats, outfile, "lat")
    NetCDF.ncwrite(lons, outfile, "lon")
    NetCDF.ncwrite(S.consts.zgeo, outfile, "height")
    NetCDF.ncwrite(S.times[1:S.tu], outfile, "time")

    # Write data variables
    for (idx, name) in enumerate(var_names)
        data_slice = data_arrays[idx][:, :, :, 1:S.tu]
        NetCDF.ncwrite(data_slice, outfile, name)
    end

    println("    netCDF file written.")
end

using NetCDF

function create_particle_track_file(cfg)
    outfile = joinpath(cfg.output_pfad, cfg.particle_track_file)
    println("Creating particle tracking file: ", outfile)
    # Remove if exists (optional)
    if isfile(outfile)
        rm(outfile)
    end
    # Create file
    ncid = NetCDF.create(outfile)
    # Define dimensions
    dim_time = NetCDF.def_dim(ncid, "time", -1)   # unlimited
    dim_particle = NetCDF.def_dim(ncid, "particle", cfg.ntrac)
    # Define coordinate variables
    var_time = NetCDF.def_var(ncid, "time", Float64, (dim_time,))
    NetCDF.put_att(ncid, var_time, "units", "seconds since 1970-01-01 00:00:00")
    var_particle = NetCDF.def_var(ncid, "particle", Int, (dim_particle,))
    NetCDF.put_att(ncid, var_particle, "long_name", "particle index")
    # Define data variables (per time, per particle)
    vars = Dict()
    vars["lon"] = NetCDF.def_var(ncid, "lon", Float32, (dim_particle, dim_time))
    NetCDF.put_att(ncid, vars["lon"], "units", "degrees_east")
    vars["lat"] = NetCDF.def_var(ncid, "lat", Float32, (dim_particle, dim_time))
    NetCDF.put_att(ncid, vars["lat"], "units", "degrees_north")
    vars["alt"] = NetCDF.def_var(ncid, "alt", Float32, (dim_particle, dim_time))
    NetCDF.put_att(ncid, vars["alt"], "units", "km")
    vars["radius"] = NetCDF.def_var(ncid, "radius", Float32, (dim_particle, dim_time))
    NetCDF.put_att(ncid, vars["radius"], "units", "nm")
    vars["is_ice"] = NetCDF.def_var(ncid, "is_ice", Int, (dim_particle, dim_time))
    NetCDF.put_att(ncid, vars["is_ice"], "long_name", "1 if ice, 0 if dust")
    vars["temp"] = NetCDF.def_var(ncid, "temp", Float32, (dim_particle, dim_time))
    NetCDF.put_att(ncid, vars["temp"], "units", "K")
    vars["pressure"] = NetCDF.def_var(ncid, "pressure", Float32, (dim_particle, dim_time))
    NetCDF.put_att(ncid, vars["pressure"], "units", "Pa")
    # Add more if needed (e.g., H2O mixing ratio at particle location)
    NetCDF.end_def(ncid)
    return ncid, var_time, vars
end

function write_particle_track_step!(ncid, var_time, vars, S, cfg, itime, current_unix_time)
    # Write time coordinate
    NetCDF.put_var(ncid, var_time, current_unix_time, start=[itime])
    # Prepare arrays for each variable
    npart = cfg.ntrac
    lon = zeros(Float32, npart)
    lat = zeros(Float32, npart)
    alt = zeros(Float32, npart)
    radius = zeros(Float32, npart)
    is_ice = zeros(Int, npart)
    temp = zeros(Float32, npart)
    pressure = zeros(Float32, npart)

    for n in 1:npart
        i = floor(Int, S.eis.xfeld[n])
        j = floor(Int, S.eis.yfeld[n])
        k = floor(Int, S.eis.zfeld[n])
        # Continuous positions
        x = S.eis.xfeld[n]
        y = S.eis.yfeld[n]
        z = S.eis.zfeld[n]
        # Convert to physical coordinates
        lon[n] = (x - 1.0f0) * 3.0f0
        lat[n] = 37.5f0 + (y - 1.0f0)
        alt[n] = 77.8f0 + (z - 1.0f0) * 0.1f0
        radius[n] = S.eis.rfeld[n]
        is_ice[n] = S.eis.rfeld[n] > S.eis.rinit[n] ? 1 : 0
        # Interpolate temperature and pressure to particle location (use nearest grid point for simplicity)
        # Since particles move continuously, we could do bilinear, but for now use nearest.
        i_int = clamp(round(Int, x), 1, cfg.igitneu)
        j_int = clamp(round(Int, y), 1, cfg.nbneu)
        k_int = clamp(round(Int, z), 1, cfg.kgitneu)
        temp[n] = S.uvwtd.tm[j_int, i_int, k_int]
        pressure[n] = S.uvwtd.dm[j_int, i_int, k_int]
    end

    # Write to netCDF (slice at current time)
    NetCDF.put_var(ncid, vars["lon"], lon, start=[1, itime], count=[npart,1])
    NetCDF.put_var(ncid, vars["lat"], lat, start=[1, itime], count=[npart,1])
    NetCDF.put_var(ncid, vars["alt"], alt, start=[1, itime], count=[npart,1])
    NetCDF.put_var(ncid, vars["radius"], radius, start=[1, itime], count=[npart,1])
    NetCDF.put_var(ncid, vars["is_ice"], is_ice, start=[1, itime], count=[npart,1])
    NetCDF.put_var(ncid, vars["temp"], temp, start=[1, itime], count=[npart,1])
    NetCDF.put_var(ncid, vars["pressure"], pressure, start=[1, itime], count=[npart,1])
    nothing
end