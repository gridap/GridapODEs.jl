module ODETools

using Gridap: GridapType
using Gridap.Helpers
using Gridap.Algebra: NonLinearSolver
using Gridap.Algebra: NonLinearOperator

export ODEOperator
import Gridap.Algebra: residual!
import Gridap.Algebra: jacobian!
import Gridap.Algebra: allocate_residual
import Gridap.Algebra: allocate_jacobian
export jacobian_unknown!
export jacobian_unknown_t!
export test_ode_operator

export ODESolver
export BackwardEuler
export test_ode_solver
import Gridap.Algebra: solve!
import Gridap.Algebra: zero_initial_guess

export ODESolution
import Gridap: solve

include("ODEOperators.jl")

include("ODESolvers.jl")

include("ODESolutions.jl")

end #module
