module ODEToolsTests

using Test

@testset "ODEOperators" begin include("ODEOperatorsTests.jl") end

@testset "ODESolvers" begin include("ODESolversTests.jl") end

end # module
