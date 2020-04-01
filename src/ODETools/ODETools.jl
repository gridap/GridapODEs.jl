module ODETools

using Gridap: GridapType
using Gridap.Helpers
using Gridap.Algebra: NonLinearSolver
using Gridap.Algebra: NonLinearOperator

export ODEOperator
export residual!
export allocate_residual
export jacobian_unknown!
export jacobian_unknown_t!
export allocate_jacobian
export test_ode_operator

export ODESolver
export BackwardEuler

include("ODEOperators.jl")

include("ODESolutions.jl")

include("ODESolvers.jl")


end #module
