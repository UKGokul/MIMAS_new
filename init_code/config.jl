# src/config.jl
# Run‑time configuration for MIMAS
# All values are taken from Fortran mimas_A_2014.f, preserving original names and types.
# This struct holds parameters that may change between runs (grid sizes, file paths, dates, switches).

struct Config
    # Grid dimensions (from #define block)
    ntrac::Int
    kgitzpr::Int
    kgitneu::Int
    kgitneu1::Int
    lunten::Int
    igitneu::Int
    igitneu1::Int
    nbneu::Int
    nbneu1::Int
    tg::Int

    # Time step parameters
    dttrans::Float64          # advection timestep (seconds)
    dttrans_2::Float64         # double timestep (used in diffusion & microphysics)

    # Run control
    nstart::Int                # 0 = cold start, 1 = restart
    i_hemi::Int                # 1 = southern hemisphere, 2 = northern
    start_month::Int
    start_day::Int
    end_month::Int
    end_day::Int
    start_year::Int             # simulation year (i_jahr)
    restart_year::String        # JAHRRESTART, used when nstart == 1
    run_info::String            # RUNINFO (log file name, optional)

    # Fixed years (from program body)
    ijahr_fixdyn::Int           # year for fixed dynamics (1976)

    # Input file paths (from #define)
    crinit::String
    h2oinit::String
    cturbkzz::String
    cbackscatter532::String
    cellipse_earth::String
    ch4_file::String
    ut1_file::String
    lima_pfad::String
    restart_pfad::String
    output_pfad::String
    stat_pfad::String
end

# Outer constructor with default values taken from the Fortran #define block and program body.
function Config(;
    ntrac = 40_000_000, # reduced for test runs
    kgitzpr = 166,
    kgitneu = 163,
    kgitneu1 = 164,
    lunten = 12,
    igitneu = 120,
    igitneu1 = 121,
    nbneu = 53,
    nbneu1 = 54,
    tg = 456,
    dttrans = 90.0,
    dttrans_2 = 180.0,
    nstart = 0,
    i_hemi = 2,
    start_month = 7,
    start_day = 10,
    end_month = 7, # is actually 8
    end_day = 13, # is actually 31
    start_year = 2008,
    restart_year = "1961",
    run_info = "run.out",
    ijahr_fixdyn = 1976,
    crinit = "/home/gu01/HPC/MIMAS/trac_3d_ini_40mio_l3.dat",
    h2oinit = "/home/gu01/HPC/MIMAS/trach2o_1dstart_zpress.dat",
    cturbkzz = "/home/gu01/HPC/knocluc3d_lima/luebkzza.dat",
    cbackscatter532 = "/home/gu01/HPC/knocluc3d_lima/backscatter532.dat",
    cellipse_earth = "/home/gu01/HPC/kmoda/ellipse_earth.dat",
    ch4_file = "/home/gu01/HPC/atmos/ch4_annual_MIMAS_0001.tab",
    ut1_file = "/home/gu01/HPC/atmos/ut1_10_may_MIMAS_MLFU.tab",
    lima_pfad = "/home/gu01/HPC/cop/CS33/LIMA-ICE/backgr/", # using CS33 instead of og one due to lack of permission
    restart_pfad = "/home/optik/dmf/CS33/LIMA-ICE/40mio/",
    output_pfad = "/home/gu01/HPC/Github/MIMAS.jl/NewTry/output/Experimental_runs/", # output file changed name as the input lima file changed, correct once after access. 
    stat_pfad = "./stat/"
)
    return Config(
        ntrac, kgitzpr, kgitneu, kgitneu1, lunten,
        igitneu, igitneu1, nbneu, nbneu1, tg,
        dttrans, dttrans_2,
        nstart, i_hemi,
        start_month, start_day, end_month, end_day,
        start_year, restart_year, run_info,
        ijahr_fixdyn,
        crinit, h2oinit, cturbkzz, cbackscatter532,
        cellipse_earth, ch4_file, ut1_file,
        lima_pfad, restart_pfad, output_pfad, stat_pfad
    )
end