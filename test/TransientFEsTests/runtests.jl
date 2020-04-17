module TransientFEToolsTests

using Test

@testset "TransientFETests" begin include("TransientFETests.jl") end

# @santiagobadia : To be eliminated when fixing the nonlinear issue
@testset "QuadraticTests" begin include("HeatEquationTests.jl") end

end # module
