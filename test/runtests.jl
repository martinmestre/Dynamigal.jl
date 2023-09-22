using GalacticDynamics
using Test, SafeTestsets
using PythonCall

include("init_python.jl")


@time begin
    @safetestset "Acceleration" begin
        include("acceleration/accelerations.jl")
    end
end