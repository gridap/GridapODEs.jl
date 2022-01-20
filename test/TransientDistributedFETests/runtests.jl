module TransientDistributedFETests

using Test

@testset "TransientDistributedCellFieldFETests" begin include("TransientDistributedCellFieldTests.jl") end

@testset "TransientMultiFieldDistributedCellFieldFETests" begin include("TransientMultiFieldDistributedCellFieldTests.jl") end

@testset "DistributedAffineFEOperatorsTests" begin include("DistributedAffineFEOperatorsTests.jl") end

@testset "DistributedBoundaryHeatEquationTests" begin include("DistributedBoundaryHeatEquationTests.jl") end

@testset "DistributedConstantFEOperatorsTests" begin include("DistributedConstantFEOperatorsTests.jl") end

@testset "DistributedDGHeatEquationTests" begin include("DistributedDGHeatEquationTests.jl") end

@testset "DistributedHeatEquationTests" begin include("DistributedHeatEquationTests.jl") end

# missing interpolate_everywhere:
#@testset "DistributedForwardEulerHeatEquationTests" begin include("DistributedForwardEulerHeatEquationTests.jl") end

# missing FESpace(triangulation,...)
#@testset "DistributedFreeSurfacePotentialFlowTests" begin include("DistributedFreeSurfacePotentialFlowTests.jl") end

# missing _jacobian(...) for DistributedDomainContribution
#@testset "DistributedHeatEquationAutoDiffTests" begin include("DistributedHeatEquationAutoDiffTests.jl") end

# missing broadcasted implementation for `A+M ≈ K`
#@testset "DistributedHeatVectorEquationTests" begin include("DistributedHeatVectorEquationTests.jl") end

end
