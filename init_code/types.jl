# src/types.jl
# =============================================================================
# Data structures for the MIMAS model state.
# (description as before)
# =============================================================================

include("constants.jl")
# table sizes from constants

# We need Config type; it will be defined in config.jl, so we assume it's available.

# ----------------------------------------------------------------------
# common_h2o
# ----------------------------------------------------------------------
mutable struct CommonH2O
    zpress_init::Vector{Float32}
    zpress_init_h2o::Vector{Float32}
    hm::Array{Float32,3}
    hminit3d::Array{Float32,3}
    hmdummy::Array{Float32,3}
    secchi::Matrix{Float32}
    hminit::Matrix{Float32}
    ellipsefak::Vector{Float32}
    oblinoa::Vector{Float32}
    true_solar_time::Vector{Float32}

    function CommonH2O(cfg::Config)
        nbneu = cfg.nbneu
        igitneu = cfg.igitneu
        kgitneu = cfg.kgitneu
        kgitzpr = cfg.kgitzpr
        new(
            zeros(Float32, kgitzpr),
            zeros(Float32, kgitzpr),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu),
            zeros(Float32, nbneu, kgitneu),
            zeros(Float32, 365),
            zeros(Float32, 365),
            zeros(Float32, 365)
        )
    end
end

# ----------------------------------------------------------------------
# common_const
# ----------------------------------------------------------------------
mutable struct CommonConst
    zgeo::Vector{Float32}
    turbkzz::Vector{Float32}
    dx::Vector{Float32}
    backradius::Vector{Float32}
    crosssec532::Vector{Float32}
    crosssec126::Vector{Float32}
    crosssec200::Vector{Float32}
    crosssec1000::Vector{Float32}
    crosssec3000::Vector{Float32}
    ext532::Vector{Float32}
    ext126::Vector{Float32}
    ext200::Vector{Float32}
    ext1000::Vector{Float32}
    ext3000::Vector{Float32}

    function CommonConst(cfg::Config)
        nbneu = cfg.nbneu
        kgitneu = cfg.kgitneu
        new(
            zeros(Float32, kgitneu),
            zeros(Float32, kgitneu),
            zeros(Float32, nbneu),
            zeros(Float32, 1500),
            zeros(Float32, 1500),
            zeros(Float32, 1500),
            zeros(Float32, 1500),
            zeros(Float32, 1500),
            zeros(Float32, 1500),
            zeros(Float32, 1500),
            zeros(Float32, 1500),
            zeros(Float32, 1500),
            zeros(Float32, 1500),
            zeros(Float32, 1500)
        )
    end
end

# ----------------------------------------------------------------------
# common_tab
# ----------------------------------------------------------------------
mutable struct CommonTab
    horiwicht_tab::Vector{Float32}
    wrichtig_tab::Vector{Float32}
    volufak_tab::Vector{Float32}
    vpsatt_tab::Vector{Float32}
    sqrt1t_tab::Vector{Float32}
    sqrtsedi_tab::Vector{Float32}
    faktkelvin_tab::Matrix{Float32}
    worka_tab::Matrix{Float32}

    function CommonTab(cfg::Config)
        nbneu = cfg.nbneu
        kgitneu = cfg.kgitneu
        new(
            zeros(Float32, nbneu),
            zeros(Float32, kgitneu),
            zeros(Float32, kgitneu),
            zeros(Float32, TAB_TMAX),
            zeros(Float32, TAB_TMAX),
            zeros(Float32, TAB_TMAX),
            zeros(Float32, TAB_RMAX, TAB_TMAX),
            zeros(Float32, TAB_RPTP, kgitneu)
        )
    end
end

# ----------------------------------------------------------------------
# common_uvwtd  (current time level fields)
# ----------------------------------------------------------------------
mutable struct CommonUVWTD
    um::Array{Float32,3}
    tm::Array{Float32,3}
    dm::Array{Float32,3}
    vm::Array{Float32,3}
    wm::Array{Float32,3}

    function CommonUVWTD(cfg::Config)
        nbneu = cfg.nbneu
        nbneu1 = cfg.nbneu1
        igitneu = cfg.igitneu
        kgitneu = cfg.kgitneu
        kgitneu1 = cfg.kgitneu1
        new(
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu1, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu1)
        )
    end
end

# ----------------------------------------------------------------------
# common_uvwtd48  (two time levels for dynamics)
# ----------------------------------------------------------------------
mutable struct CommonUVWTD48
    um1::Array{Float32,3}
    um2::Array{Float32,3}
    tm1::Array{Float32,3}
    tm2::Array{Float32,3}
    dm1::Array{Float32,3}
    dm2::Array{Float32,3}
    zm1::Array{Float32,3}
    zm2::Array{Float32,3}
    vm1::Array{Float32,3}
    vm2::Array{Float32,3}
    wm1::Array{Float32,3}
    wm2::Array{Float32,3}

    function CommonUVWTD48(cfg::Config)
        nbneu = cfg.nbneu
        nbneu1 = cfg.nbneu1
        igitneu = cfg.igitneu
        kgitneu = cfg.kgitneu
        kgitneu1 = cfg.kgitneu1
        new(
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu),
            zeros(Float32, nbneu1, igitneu, kgitneu),
            zeros(Float32, nbneu1, igitneu, kgitneu),
            zeros(Float32, nbneu, igitneu, kgitneu1),
            zeros(Float32, nbneu, igitneu, kgitneu1)
        )
    end
end

# ----------------------------------------------------------------------
# common_eis  (particle data)
# ----------------------------------------------------------------------
mutable struct CommonEIS
    xfeld::Vector{Float32}
    yfeld::Vector{Float32}
    zfeld::Vector{Float32}
    rfeld::Vector{Float32}
    rinit::Vector{Float32}

    function CommonEIS(cfg::Config)
        ntrac = cfg.ntrac
        new(
            zeros(Float32, ntrac),
            zeros(Float32, ntrac),
            zeros(Float32, ntrac),
            zeros(Float32, ntrac),
            zeros(Float32, ntrac)
        )
    end
end

# ----------------------------------------------------------------------
# common_netcdf  (output buffers)
# ----------------------------------------------------------------------
mutable struct CommonNetCDF
    aq::Array{Float32,4}
    bq::Array{Float32,4}
    cq::Array{Float32,4}
    dq::Array{Float32,4}
    eq::Array{Float32,4}
    fq::Array{Float32,4}
    gq::Array{Float32,4}
    hq::Array{Float32,4}
    iq::Array{Float32,4}
    jq::Array{Float32,4}
    kq::Array{Float32,4}
    lq::Array{Float32,4}
    mq::Array{Float32,4}
    nq::Array{Float32,4}
    oq::Array{Float32,4}
    pq::Array{Float32,4}

    function CommonNetCDF(cfg::Config)
        kgitneu = cfg.kgitneu
        igitneu = cfg.igitneu
        nbneu = cfg.nbneu
        tg = cfg.tg
        dims = (kgitneu, igitneu, nbneu, tg)
        new(
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims),
            zeros(Float32, dims)
        )
    end
end

# ----------------------------------------------------------------------
# Top‑level MIMASState – contains all common blocks and additional work arrays
# ----------------------------------------------------------------------
mutable struct MIMASState
    h2o::CommonH2O
    consts::CommonConst
    tab::CommonTab
    uvwtd::CommonUVWTD          # current interpolated fields
    uvwtd48::CommonUVWTD48       # two time levels
    eis::CommonEIS
    netcdf::CommonNetCDF

    # Derived arrays (from Fortran main program)
    j_mesopause::Vector{Int}     # (nbneu)
    j_zpress86km::Vector{Int}    # (nbneu)
    j_zbound_ob::Vector{Int}     # (nbneu)  – upper boundary for particle rebirth

    # Work arrays (to avoid allocations in time loops)
    workw::Array{Float32,3}       # same as wm (nbneu, igitneu, kgitneu1)
    change_n::Vector{Float32}     # (ntrac)  – per‑particle H2O change
    work1::Array{Float32,3}       # (nbneu, igitneu, kgitneu) – temporary for microphysics
    change::Array{Float32,3}      # (nbneu, igitneu, kgitneu) – accumulated H2O change

    # Pre‑computed factors (optional, for performance)
    ax::Vector{Float32}            # dttrans_2 / dx  (nbneu)

    # Particle‑distribution helpers
    bbreite::Vector{Float32}       # (35) latitude boundaries
    bsum::Vector{Float32}          # (35) cumulative sums
    sfactor::Float32                # sin(55°)
    invfactor::Float32              # 1 - sfactor

    # Output counters
    tu::Int                         # number of outputs written
    times::Vector{Float32}          # unix times of outputs (length tg)

    function MIMASState(cfg::Config)
        h2o = CommonH2O(cfg)
        consts = CommonConst(cfg)
        tab = CommonTab(cfg)
        uvwtd = CommonUVWTD(cfg)
        uvwtd48 = CommonUVWTD48(cfg)
        eis = CommonEIS(cfg)
        netcdf = CommonNetCDF(cfg)

        nbneu = cfg.nbneu
        igitneu = cfg.igitneu
        kgitneu = cfg.kgitneu
        kgitneu1 = cfg.kgitneu1
        ntrac = cfg.ntrac

        j_mesopause = zeros(Int, nbneu)
        j_zpress86km = zeros(Int, nbneu)
        j_zbound_ob = zeros(Int, nbneu)

        workw = zeros(Float32, nbneu, igitneu, kgitneu1)
        change_n = zeros(Float32, ntrac)
        work1 = zeros(Float32, nbneu, igitneu, kgitneu)
        change = zeros(Float32, nbneu, igitneu, kgitneu)

        bbreite = zeros(Float32, 35)
        bsum = zeros(Float32, 35)
        sfactor = 0.0f0
        invfactor = 0.0f0

        ax = zeros(Float32, nbneu)

        tu = 0
        times = zeros(Float32, cfg.tg)

        return new(
            h2o, consts, tab, uvwtd, uvwtd48, eis, netcdf,
            j_mesopause, j_zpress86km, j_zbound_ob,
            workw, change_n, work1, change,
            ax, bbreite, bsum, sfactor, invfactor,
            tu, times
        )
    end
end

# ----------------------------------------------------------------------
# Optional: helper function to update ax after dx is known
function update_ax!(S::MIMASState, cfg::Config)
    @. S.ax = cfg.dttrans_2 / S.consts.dx
end
# ----------------------------------------------------------------------