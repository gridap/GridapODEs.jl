module ODETools

using Test

using Gridap: GridapType
using Gridap.Helpers
using Gridap.Algebra: NonlinearSolver
using Gridap.Algebra: LinearSolver
using Gridap.Algebra: NonlinearOperator
using Gridap.Algebra: AffineOperator

export ODEOperator
using Gridap.Algebra: residual
using Gridap.Algebra: jacobian
using Gridap.Algebra: symbolic_setup
using Gridap.Algebra: numerical_setup
using Gridap.Algebra: numerical_setup!
import Gridap.Algebra: residual!
import Gridap.Algebra: jacobian!
import Gridap.Algebra: allocate_residual
import Gridap.Algebra: allocate_jacobian
export allocate_cache
export update_cache!
export jacobian!
export jacobian_t!
export test_ode_operator

export ODESolver
export solve_step!
export test_ode_solver
import Gridap: solve
import Gridap.Algebra: solve!
import Gridap.Algebra: zero_initial_guess

export BackwardEuler
export ForwardEuler
export MidPoint
export ThetaMethod
import Gridap.Algebra: fill_entries!

export ODESolution
export test_ode_solution
import Base: iterate

include("ODEOperators.jl")

include("ODESolvers.jl")

include("ODESolutions.jl")

end #module
