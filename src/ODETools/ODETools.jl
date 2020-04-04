module ODETools

using Test

using Gridap: GridapType
using Gridap.Helpers
using Gridap.Algebra: NonlinearSolver
using Gridap.Algebra: NonlinearOperator

export ODEOperator
import Gridap.Algebra: residual!
import Gridap.Algebra: jacobian!
import Gridap.Algebra: allocate_residual
import Gridap.Algebra: allocate_jacobian
export jacobian!
export jacobian_t!
export test_ode_operator

export ODESolver
export BackwardEuler
export test_ode_solver
import Gridap.Algebra: solve!
import Gridap.Algebra: zero_initial_guess

export ODESolution
import Gridap: solve
export test_ode_solution

include("ODEOperators.jl")

include("ODESolvers.jl")

include("ODESolutions.jl")

end #module
