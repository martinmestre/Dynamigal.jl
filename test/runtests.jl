using GalacticDynamics
using Test, SafeTestsets

@time begin
    @safetestset "Acceleration" begin
        include("acceleration/acceleration.jl")
    end
end