module TransientFEToolsTests

using Test

@testset "TransientFETests" begin include("TransientFETests.jl") end

@testset "TransientFEOperatorsTests" begin include("TransientFEOperatorsTests.jl") end

@testset "HeatEquationTests" begin include("HeatEquationTests.jl") end

@testset "HeatVectorEquationTests" begin include("HeatVectorEquationTests.jl") end

@testset "VectorHeatEquationTests" begin include("VectorHeatEquationTests.jl") end

@testset "StokesEquationTests" begin include("StokesEquationTests.jl") end

end # module
