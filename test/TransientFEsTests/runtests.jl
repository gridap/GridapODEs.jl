module TransientFEToolsTests

using Test

@testset "TransientFETests" begin include("TransientFETests.jl") end

@testset "HeatEquationTests" begin include("HeatEquationTests.jl") end

end # module
