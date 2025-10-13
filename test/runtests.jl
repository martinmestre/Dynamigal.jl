using Dynamigal
using Test
using PythonCall
using BenchmarkTools


pyimport("sys")."path".append("")
pyimport("sys")."path".append("python")
accelerations_py = pyimport("accelerations")
au = pyimport("astropy.units")
gd = pyimport("gala.dynamics")
gp = pyimport("gala.potential")
gu = pyimport("gala.units")
gi = pyimport("gala.integrate")

# @time begin
#     include("acceleration/test_accelerations.jl")
# end
# @time begin
#     include("ode/test_odes.jl")
# end
@time begin
    include("orbit/test_orbits.jl")
end

