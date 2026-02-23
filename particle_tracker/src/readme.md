# MIMAS.jl

Julia implementation of the MIMAS (Mesospheric Ice Microphysics and Advection Scheme) model for simulating noctilucent clouds (NLCs), including particle microphysics, transport, photochemistry, and radiative properties.

This version reproduces the structure and functionality of the original Fortran model while adopting a modular, maintainable Julia architecture.

---

# Project Structure

```
MIMAS.jl/
├── src/
│   ├── constants.jl
│   ├── config.jl
│   ├── types.jl
│   ├── init.jl
│   ├── io.jl
│   ├── subroutines.jl
│   └── main.jl
├── test/
└── README.md
```

---

# File Descriptions

## constants.jl

Contains:

* Fundamental physical constants (π, Earth radius, Boltzmann constant, etc.)
* Fixed grid dimensions
* Parameters that never change between runs

These values define the physical framework of the model and should generally not be modified.

---

## config.jl

Contains all user-adjustable simulation settings.

This is the only file that typical users need to edit before running a simulation.

It defines:

* Grid size
* Time stepping
* Simulation period
* Hemisphere selection
* Input/output file paths
* Run switches

Full details are provided in the Config Section below.

---

## types.jl

Defines the model state structure:

```
mutable struct MIMASState
```

This replaces the original Fortran COMMON blocks.

All particle arrays, Eulerian fields, and work arrays are allocated and initialised to zero when the state is created.

---

## init.jl

Initialisation routines that:

* Compute grid geometry (dx, zgeo)
* Compute auxiliary indices:

  * j_mesopause
  * j_zpress86km
  * j_zbound_ob
* Build lookup tables:

  * Saturation vapour pressure
  * Kelvin factor
  * Other thermodynamic tables

This stage prepares the model state before entering the time loop.

---

## io.jl

Handles all file I/O operations.

### Low-level routines

* read_fortran_record
* write_fortran_record

Used for reading/writing Fortran unformatted binary files.

### High-level routines

* Read global initialisation files
* Read dynamics snapshots
* Write restart files
* Write netCDF output

---

## subroutines.jl

Contains all physics modules:

* Advection
* Diffusion
* Microphysics
* Photolysis
* Zenith angle calculation
* Particle birth control
* Utility functions

This file contains the scientific core of the model.

---

## main.jl

The main simulation driver.

Execution flow:

1. Load configuration
2. Create model state
3. Run initialisation
4. Enter time loop:

   * Months
   * Days
   * Hours
   * 20 inner timesteps (90 s)

Output is written every 6 hours.

---

# Config Struct (Detailed Explanation)

All user-modifiable parameters are defined in:

```
struct Config
```

---

## Grid Dimensions

| Field    | Description                            | Default                 |
| -------- | -------------------------------------- | ----------------------- |
| ntrac    | Number of particles                    | 40,000,000 (production) |
| kgitzpr  | Pressure levels for H2O initialisation | 166                     |
| kgitneu  | Geometric height levels                | 163                     |
| kgitneu1 | kgitneu + 1 (vertical wind)            | 164                     |
| lunten   | Lowest vertical transport index        | 12                      |
| igitneu  | Longitude points                       | 120                     |
| igitneu1 | igitneu + 1                            | 121                     |
| nbneu    | Latitude points                        | 53                      |
| nbneu1   | nbneu + 1                              | 54                      |
| tg       | Number of output times                 | 456                     |

If grid dimensions are modified, all input files must match the new grid.

---

## Time Step Parameters

| Field     | Description        | Default |
| --------- | ------------------ | ------- |
| dttrans   | Advection timestep | 90 s    |
| dttrans_2 | Double timestep    | 180 s   |

---

## Run Control

| Field                  | Description                           |
| ---------------------- | ------------------------------------- |
| nstart                 | 0 = cold start, 1 = restart           |
| i_hemi                 | 1 = Southern Hemisphere, 2 = Northern |
| start_month, start_day | Simulation start date                 |
| end_month, end_day     | Simulation end date                   |
| start_year             | Simulation year                       |
| restart_year           | Restart file year                     |
| run_info               | Log filename                          |

---

## Fixed Year Parameter

| Field        | Description                    |
| ------------ | ------------------------------ |
| ijahr_fixdyn | Year for fixed dynamics (1976) |

---

## Input File Paths

All paths must exist.

| Field           | Description                    |
| --------------- | ------------------------------ |
| crinit          | Restart/initial particle file  |
| h2oinit         | H2O pressure profile           |
| cturbkzz        | Turbulence coefficient file    |
| cbackscatter532 | Backscatter cross-section file |
| cellipse_earth  | Earth orbit parameters         |
| ch4_file        | Methane time series            |
| ut1_file        | Unix time reference            |
| lima_pfad       | Dynamics snapshot directory    |
| restart_pfad    | Restart directory              |
| output_pfad     | NetCDF output directory        |
| stat_pfad       | Statistics directory           |

---

# Running a Simulation

## Step 1: Edit Config

Modify:

```
src/config.jl
```

Set particle count, dates, and file paths.

---

## Step 2: Run

From project root:

```
julia --project=. src/main.jl
```

Monitor terminal output for Courant numbers and timestep diagnostics.

---

## Output

After completion:

```
output_pfad/<start_year>.nc
```

Example:

```
2008.nc
```

---

# Example: 1-Day Test Run

In config.jl:

```julia
ntrac = 1000
start_day = 10
end_day = 10
```

Run:

```
julia --project=. src/main.jl
```

Should finish in a few minutes.

---

# Output NetCDF Variables

Dimensions: (height, lon, lat, time)

| Variable  | Units          | Description                  |
| --------- | -------------- | ---------------------------- |
| beta532   | 10⁻¹⁰ m⁻¹ sr⁻¹ | Volume backscatter at 532 nm |
| densice   | g km⁻³         | Ice water density            |
| r_ice     | nm             | Mean ice particle radius     |
| n_ice     | cm⁻³           | Ice number density           |
| n_dust    | cm⁻³           | Dust number density          |
| r_ice_eff | nm             | Effective radius             |
| H2O       | ppmv           | Water vapour mixing ratio    |
| zon_win   | m/s            | Zonal wind                   |
| mer_win   | m/s            | Meridional wind              |
| tem       | K              | Temperature                  |
| press_alt | km             | Pressure altitude            |
| ext_126   | m⁻¹            | Extinction at 126 nm         |
| ext_200   | m⁻¹            | Extinction at 200 nm         |
| ext_532   | m⁻¹            | Extinction at 532 nm         |
| ext_1000  | m⁻¹            | Extinction at 1000 nm        |
| ext_3000  | m⁻¹            | Extinction at 3000 nm        |

---

## Coordinates

| Coordinate | Description                   |
| ---------- | ----------------------------- |
| lat        | 37.5 … 37.5+(nbneu−1) (°N)    |
| lon        | 0, 3, 6, … , 357 (°E)         |
| height     | 77.8 km upward (Δz = 0.1 km)  |
| time       | Unix seconds since 1970‑01‑01 |

---

# Performance Tips

## Production Runs (40M particles)

* Serial execution will be slow.
* Add threading to particle loops:

```julia
Threads.@threads for ...
```

## Optimisation Suggestion

* Preallocate temporary arrays in sub_transwalcek!
* Avoid allocations inside loops

## I/O

* Store input files on a fast local disk
* Avoid network-mounted drives

---

# Restart Capability (Experimental)

If:

```
nstart = 0
```

A restart file is written after initialisation.

If:

```
nstart = 1
```

The restart file is read.

Restart functionality is not fully validated.

---

# Troubleshooting

## File Not Found

Double-check all paths in config.jl.

---

## Record Size Mismatch

Binary files are assumed to:

* Be Fortran unformatted
* Use big-endian
* Have 4-byte record headers

Modify read_fortran_record if your files differ.

---

## NetCDF Errors

Ensure NetCDF.jl is installed:

```
] add NetCDF
```

If you encounter:

```
UndefVarError: def_dim
```

Your NetCDF.jl API version may differ. Update write_netcdf! accordingly.

---

# Notes for Developers

* The model structure mirrors the original Fortran code.
* COMMON blocks are replaced with structured state objects.
* Particle transport is Lagrangian.
* Water vapour transport is Eulerian.
* Both use the same eddy diffusivity (Kz).

---

# License

Specify project license here.

---

# Author

Specify author information here.
