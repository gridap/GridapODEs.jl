module ODETools

using Gridap: GridapType
using Gridap.Helpers
using Gridap.Algebra: NonLinearSolver
using Gridap.Algebra: NonLinearOperator

import Gridap: solve
import Gridap: solve!

export ODEOperator
export residual!
export allocate_residual
export jacobian_unknown!
export jacobian_unknown_t!
export allocate_jacobian
export test_ode_operator

export ODESolver
export BackwardEuler
export test_ode_solver


include("ODEOperators.jl")

include("ODESolvers.jl")

include("ODESolutions.jl")



end #module
