module TransientDistributedFETests

using Test

@testset "TransientDistributedCellFieldFETests" begin include("TransientDistributedCellFieldTests.jl") end

@testset "DistributedAffineFEOperatorsTests" begin include("DistributedAffineFEOperatorsTests.jl") end

@testset "DistributedBoundaryHeatEquationTests" begin include("DistributedBoundaryHeatEquationTests.jl") end

@testset "DistributedConstantFEOperatorsTests" begin include("DistributedConstantFEOperatorsTests.jl") end

end
