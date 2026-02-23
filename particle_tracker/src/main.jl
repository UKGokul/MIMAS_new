include("constants.jl")
include("config.jl")
include("types.jl")
include("init.jl")
include("subroutines.jl")
include("io.jl")

using Random
using NetCDF

println("BOOT OK")
println("="^50)

# ----------------------------------------------------------------------
# 1. Configuration
# ----------------------------------------------------------------------
cfg = Config()

println("Configuration loaded:")
println("  Grid: nbneu=$(cfg.nbneu), igitneu=$(cfg.igitneu), kgitneu=$(cfg.kgitneu)")
println("  Start: $(cfg.start_year)-$(cfg.start_month)-$(cfg.start_day)")
println("  End: $(cfg.end_month)-$(cfg.end_day)")
println("  Hemisphere: $(cfg.i_hemi == 1 ? "Southern" : "Northern")")
println("  Start type: $(cfg.nstart == 0 ? "Cold start" : "Restart")")
println("  Output path: $(cfg.output_pfad)")
println("="^50)

# Ensure output directory exists
mkpath(cfg.output_pfad)

# ----------------------------------------------------------------------
# 2. Main simulation function
# ----------------------------------------------------------------------
function run_simulation(cfg)
    println("\n=== Starting simulation initialisation ===")

    # ------------------------------------------------------------------
    # 2.1 Allocate state and geometry
    # ------------------------------------------------------------------
    S = allocate_state(cfg)
    init_constants!(S, cfg)

    # ------------------------------------------------------------------
    # 2.2 Read global initialisation files
    # ------------------------------------------------------------------
    xlyman_obs, scale_ch4, unix_time = init_global!(S, cfg, cfg.start_year)

    # Quick size check
    println("Size of hm: ", size(S.h2o.hm))
    println("Size of um1: ", size(S.uvwtd48.um1))
    println("Size of xfeld: ", length(S.eis.xfeld))
    println("Size of aq: ", size(S.netcdf.aq))
    println("Types OK")

    # ------------------------------------------------------------------
    # 2.3 Month lengths for start year
    # ------------------------------------------------------------------
    ntage_monat, max_doy = set_monate(cfg.start_year)
    println("Days per month: ", ntage_monat)
    println("Max day of year: ", max_doy)

    # ------------------------------------------------------------------
    # 2.4 Read H2O initialisation
    # ------------------------------------------------------------------
    println("\n--- Reading H2O initialisation ---")
    read_h2o_init!(S, cfg)

    # ------------------------------------------------------------------
    # 2.5 Determine start time for first dynamics snapshot
    # ------------------------------------------------------------------
    println("\n--- Determining start time for dynamics ---")
    i_jahrst = cfg.start_year
    i_monatst = cfg.start_month
    i_tagst = cfg.start_day - 1
    i_stdst = 24
    if cfg.start_day == 1
        i_monatst -= 1
        if i_monatst == 0
            i_monatst = 12
            i_jahrst -= 1
        end
        ntage_monat_temp, _ = set_monate(i_jahrst)
        i_tagst = ntage_monat_temp[i_monatst]
    end
    println("  First dynamics snapshot: $(i_jahrst)-$(i_monatst)-$(i_tagst) $(i_stdst):00")

    # ------------------------------------------------------------------
    # 2.6 Read first dynamics snapshot into *1
    # ------------------------------------------------------------------
    println("\n--- Reading first dynamics snapshot ---")
    read_dynamics!(S, cfg, cfg.i_hemi, cfg.ijahr_fixdyn, i_jahrst, i_monatst, i_tagst, i_stdst)

    # ------------------------------------------------------------------
    # 2.7 Check Courant numbers
    # ------------------------------------------------------------------
    println("\n--- Checking Courant numbers ---")
    sub_maxmin!(S.uvwtd48.um1, S.uvwtd48.vm1, S.uvwtd48.wm1, S, cfg)

    # ------------------------------------------------------------------
    # 2.8 Compute mesopause and 86 km indices
    # ------------------------------------------------------------------
    println("\n--- Computing mesopause and 86 km indices ---")
    compute_mesopause_and_86km!(S, cfg)
    compute_j_zbound_ob!(S, cfg)

    # ------------------------------------------------------------------
    # 2.9 Initialise particles
    # ------------------------------------------------------------------
    println("\n--- Initialising particles ---")
    init_staub!(S, cfg)

    # ------------------------------------------------------------------
    # 2.10 Write restart file (if cold start)
    # ------------------------------------------------------------------
    if cfg.nstart == 0
        write_restart!(S, cfg)
    end

    # ------------------------------------------------------------------
    # 2.11 Compute lookup tables
    # ------------------------------------------------------------------
    println("\n--- Computing lookup tables ---")
    init_tables!(S, cfg)

    # Store initial Unix time for output
    current_unix_time = unix_time

    # ------------------------------------------------------------------
    # 3. Main time loop (months, days, hours, inner steps)
    # ------------------------------------------------------------------
    println("\n" * "="^50)
    println("=== Starting main time loop ===")

    # Compute initial day-of-year for the day before start
    itag_im_jahr = 0
    for im in 1:cfg.start_month-1
        itag_im_jahr += ntage_monat[im]
    end
    itag_im_jahr += cfg.start_day
    itag_im_jahr -= 1   # day before start
    println("Initial itag_im_jahr (day before start): ", itag_im_jahr)
    first_timestep = true

    # Month loop
    i_monat2 = cfg.start_month
    i_jahr = cfg.start_year
    while i_monat2 <= cfg.end_month
        i_monat = i_monat2
        if i_monat2 > 12
            i_monat = i_monat2 - 12
            i_jahr += 1
        end

        ntage_monat_local, max_doy_local = set_monate(i_jahr)

        monanf = 1
        monend = ntage_monat_local[i_monat]
        if i_monat2 == cfg.start_month
            monanf = cfg.start_day
        end
        if i_monat2 == cfg.end_month
            monend = cfg.end_day
        end

        println("  Month: ", i_monat, " (", i_monat2, "), year: ", i_jahr,
                ", days: ", monanf, " to ", monend)

        # Day loop
        for i_tag in monanf:monend
            itag_im_jahr += 1
            if itag_im_jahr > max_doy_local
                itag_im_jahr = 1
            end
            ntage_monat_local, max_doy_local = set_monate(i_jahr)

            println("    itag_im_jahr: ", itag_im_jahr,
                    " | date: ", i_jahr, "-", i_monat, "-", i_tag)

            # --- Fortran: srand(itag_im_jahr), reset ntg, hour loop ---
            Random.seed!(itag_im_jahr)       # seed the random generator for the day
            ntg = 0                           # inner timestep counter

            # Hour loop (1 to 24)
            for i_std in 1:24
                # Copy *1 fields to working arrays (start of hour)
                S.uvwtd.um .= S.uvwtd48.um1
                S.uvwtd.tm .= S.uvwtd48.tm1
                S.uvwtd.dm .= S.uvwtd48.dm1
                S.uvwtd.vm .= S.uvwtd48.vm1
                S.uvwtd.wm .= S.uvwtd48.wm1

                # Compute hminit3d from current pressure altitude (zm1)
                h2oinit_zpr_to_zgeo!(S, cfg)

                # Apply scaling factors (0.9 and CH4 scaling)
                S.h2o.hminit3d .*= 0.9f0
                S.h2o.hminit3d .*= scale_ch4

                # If this is the first timestep of the whole simulation, set hm = hminit3d
                if first_timestep
                    S.h2o.hm .= S.h2o.hminit3d
                    first_timestep = false
                end

                # --- Fortran block for initial photolysis adjustment ---
                # Save current hm
                S.h2o.hmdummy .= S.h2o.hm
                # Set hm to hminit3d * 1.2
                S.h2o.hm .= S.h2o.hminit3d .* 1.2f0
                # Call sub_zenit! with ntg=1
                sub_zenit!(1, itag_im_jahr, cfg.i_hemi, S, cfg)
                # Call photolysis with mnudge=0, large timestep, constant Lyman‑α
                xlyman_const_arr = fill(4.46597f0, 366)   # constant array for this call
                sub_photolyse!(0, 100000.0f0, itag_im_jahr, cfg.i_hemi,
                               xlyman_const_arr, S, cfg)
                # Update hminit3d with the result and restore hm
                S.h2o.hminit3d .= S.h2o.hm
                S.h2o.hm .= S.h2o.hmdummy

                # (Optional: compute zonal mean hminit, as in Fortran)
                for k in 1:cfg.kgitneu
                    for j in 1:cfg.nbneu
                        dsum = 0.0f0
                        for i in 1:cfg.igitneu
                            dsum += S.h2o.hminit3d[j,i,k]
                        end
                        S.h2o.hminit[j,k] = dsum / cfg.igitneu
                    end
                end

                # Read next dynamics snapshot into *2 for this hour
                read_dynamics_into!(S.uvwtd48.um2, S.uvwtd48.tm2, S.uvwtd48.dm2,
                                    S.uvwtd48.vm2, S.uvwtd48.wm2, S.uvwtd48.zm2,
                                    S, cfg, cfg.i_hemi, cfg.ijahr_fixdyn,
                                    i_jahr, i_monat, i_tag, i_std)

                # Smooth wm2
                smooth_wm!(S.uvwtd48.wm2, cfg)

                # Check Courant numbers for *2
                sub_maxmin!(S.uvwtd48.um2, S.uvwtd48.vm2, S.uvwtd48.wm2, S, cfg)

                # Inner timestep loop (20 iterations)
                for idummy in 1:20
                    ntg += 1
                    # Interpolation weights (a1 for *1, a2 for *2)
                    a1 = 1.0f0 - (ntg % 20) / 20.0f0
                    a2 = 1.0f0 - a1

                    # Interpolate working arrays
                    S.uvwtd.um .= a1 .* S.uvwtd48.um1 .+ a2 .* S.uvwtd48.um2
                    S.uvwtd.tm .= a1 .* S.uvwtd48.tm1 .+ a2 .* S.uvwtd48.tm2
                    S.uvwtd.dm .= a1 .* S.uvwtd48.dm1 .+ a2 .* S.uvwtd48.dm2
                    S.uvwtd.vm .= a1 .* S.uvwtd48.vm1 .+ a2 .* S.uvwtd48.vm2
                    S.uvwtd.wm .= a1 .* S.uvwtd48.wm1 .+ a2 .* S.uvwtd48.wm2

                    # First transport step
                    S.h2o.hm .= max.(0.01f0, S.h2o.hm)
                    sub_transwalcek!(S, cfg)

                    # Second transport step
                    S.h2o.hm .= max.(0.01f0, S.h2o.hm)
                    sub_transwalcek!(S, cfg)

                    # Diffusion
                    sub_diffu_h2o!(S, cfg)

                    # Zenith angle (updates secchi)
                    sub_zenit!(ntg, itag_im_jahr, cfg.i_hemi, S, cfg)

                    # Photolysis (mnudge=1, delta_t = dttrans_2)
                    sub_photolyse!(1, cfg.dttrans_2, itag_im_jahr,
                                   cfg.i_hemi, xlyman_obs, S, cfg)

                    # Particle microphysics and transport
                    sub_tracertransp!(S, cfg)
                    update_mesopause!(S, cfg)
                    particle_birth_control!(S, cfg, ntg)
                end

                # After inner loop, swap *2 into *1 for next hour
                S.uvwtd48.um1 .= S.uvwtd48.um2
                S.uvwtd48.tm1 .= S.uvwtd48.tm2
                S.uvwtd48.dm1 .= S.uvwtd48.dm2
                S.uvwtd48.vm1 .= S.uvwtd48.vm2
                S.uvwtd48.wm1 .= S.uvwtd48.wm2
                S.uvwtd48.zm1 .= S.uvwtd48.zm2

                # Output every 6 hours
                if i_std in (6, 12, 18, 24)
                    schreibe_beta0!(S, cfg, current_unix_time)
                    current_unix_time += 21600.0f0
                end
            end # hour loop
        end # day loop
        i_monat2 += 1
    end # month loop

    println("=== End of simulation ===")
    # Write netCDF output
    write_netcdf!(S, cfg)
    println("="^50)
end

# ----------------------------------------------------------------------
# 4. Run the simulation
# ----------------------------------------------------------------------
run_simulation(cfg)