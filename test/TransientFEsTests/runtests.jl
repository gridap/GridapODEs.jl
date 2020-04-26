module TransientFEToolsTests

using Test

@testset "TransientFETests" begin include("TransientFETests.jl") end

@testset "HeatEquationTests" begin include("HeatEquationTests.jl") end

@testset "HeatVectorEquationTests" begin include("HeatVectorEquationTests.jl") end

end # module
