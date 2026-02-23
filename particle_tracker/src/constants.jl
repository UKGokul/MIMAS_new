# src/constants.jl
# Physical and numerical constants for MIMAS
# All values are taken from Fortran mimas_A_2014.f, preserving original names and types.
# Units are indicated where known.

# --------------------------------------------
# Grid dimensions (from #define)
# --------------------------------------------

const ntrac   = 40_000_000          # number of particles (40 million)
const kgitzpr = 166                 # number of pressure levels for H2O initialisation
const kgitneu = 163                 # number of vertical levels (geometric height)
const kgitneu1 = 164                # kgitneu + 1 (for vertical wind)
const lunten  = 12                  # lowest vertical index for transport (starting level)
const igitneu = 120                 # number of longitude points
const igitneu1 = 121                # igitneu + 1 (used in advection)
const nbneu   = 53                  # number of latitude points
const nbneu1  = 54                  # nbneu + 1 (for meridional wind)
const tg      = 456                 # number of output time steps (for netCDF)

# --------------------------------------------
# Physical constants
# --------------------------------------------

const pi       = 3.14159265358979f0   # π
const dttrans  = 90.0f0               # advection timestep (seconds) – 960 steps per day
const dttrans_2 = 180.0f0             # double timestep (used in diffusion and microphysics)
const rearth   = 6_366_197.0f0        # Earth radius (m)
const gmeso    = 9.55f0               # gravity in mesosphere (m/s²) – approximate constant

const xkbolz   = 1.3805f-23           # Boltzmann constant (J/K)
const xmh2o    = 2.99f-26             # mass of H2O molecule (kg)
const xmair    = 4.845f-26            # mean mass of air molecule (kg)
const rhoice   = 932.0f0              # density of ice (kg/m³)

# --------------------------------------------
# Numerical limits and table sizes
# --------------------------------------------

const NMAX     = 500                  # maximum size for tridiagonal solver arrays
const TAB_RPTP = 200                  # number of radius bins for particle temperature table (1–200 nm)
const TAB_RMAX = 300                  # number of radius bins for Kelvin factor table (0.1–30 nm in 0.1 steps)
const TAB_TMAX = 71                   # number of temperature bins for tables (100–170 K in 1 K steps)
const beta05   = 0.24f0               # empirical factor for particle temperature (12.0/50.0)
const t82km    = 13.0f0               # reference temperature at 82 km (K) – used in worka_tab

# --------------------------------------------
# Geometric constants (from common_const PARAMETER)
# --------------------------------------------

const dz = 100.0f0                    # vertical grid spacing (m) – actually 0.1 km, but stored as 100.0 m
const dy = 111_111.1111f0             # meridional grid spacing (m) – 1° latitude in metres

# --------------------------------------------
# Additional constants found in code (used in multiple subroutines)
# --------------------------------------------

# Chapman function coefficients (sub_zenit)
const skh    = 7.0f0                  # scale height factor
const zgeo_86 = 86.0f0                # reference height (km) for zenith angle calculation
const za     = 1.0606963f0
const zb     = 0.55643831f0
const zd     = 1.7245609f0
const zccc   = 1.0619896f0
const zf     = 0.56498823f0
const zg     = 0.06651874f0

# Photodissociation cross section for H2O (cm²) – from sub_photolyse
const photodiss_cross_h2o = 1.53f-17

# Fractional release factor for methane – from sub_init_global
const fr_alpha = 0.95f0

# Small offset for boundary handling – from sub_tracertransp
const eps = 1.0f-4

