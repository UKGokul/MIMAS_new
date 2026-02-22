# test_netcdf.jl
include("constants.jl")
include("config.jl")
include("types.jl")
include("io.jl")   # this contains write_netcdf!

using NetCDF

# Create a dummy configuration (use your actual paths, but output path can be temporary)
cfg = Config(; output_pfad = "./test_output")  # use a local test directory
mkpath(cfg.output_pfad)

# Create a dummy state S with minimal fields needed for write_netcdf!
# We need S.netcdf with 16 4D arrays, S.consts.zgeo, S.tu, S.times.
S = MIMASState(cfg)   # this allocates all arrays (zeros)

# Fill a few dummy values to make the output interesting
S.consts.zgeo .= 77.8f0 .+ (0:cfg.kgitneu-1)*0.1f0
S.tu = 2   # pretend we have 2 output times
S.times[1] = Float32(1.21e9)
S.times[2] = Float32(1.21e9 + 21600)
# Fill some data arrays with simple patterns
for k in 1:cfg.kgitneu
    for i in 1:cfg.igitneu
        for j in 1:cfg.nbneu
            S.netcdf.aq[k,i,j,1] = 0.1f0 * k
            S.netcdf.aq[k,i,j,2] = 0.2f0 * k
            S.netcdf.bq[k,i,j,1] = 0.3f0 * i
            S.netcdf.bq[k,i,j,2] = 0.4f0 * i
            # ... fill others if you want, but not necessary for testing
        end
    end
end

# Now try to write the netCDF file
println("Testing netCDF write...")
write_netcdf!(S, cfg)
println("File written. Check ./test_output/")