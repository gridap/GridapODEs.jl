module GridapTimeStepperRunTests

using Test

@time @testset "ODETools" begin include("ODEsTests/runtests.jl") end

# include("../bench/runbenchs.jl")

end #module
