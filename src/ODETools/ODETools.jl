module ODETools

using Gridap: GridapType
using Gridap.Helpers
using Gridap.Algebra: NonLinearSolver
using Gridap.Algebra: NonLinearOperator

import Gridap: solve
import Gridap: solve!
using Gridap.Algebra: residual
import Gridap.Algebra: residual!
import Gridap.Algebra: jacobian!
import Gridap.Algebra: allocate_residual
import Gridap.Algebra: allocate_jacobian
import Gridap.Algebra: zero_initial_guess
using Gridap.Algebra: jacobian

export ODEOperator
# export residual!
# export allocate_residual
export jacobian_unknown!
export jacobian_unknown_t!
export allocate_jacobian
export test_ode_operator

export ODESolver
export BackwardEuler
export test_ode_solver
import Gridap.Algebra: solve!

export ODESolution


include("ODEOperators.jl")

include("ODESolvers.jl")

include("ODESolutions.jl")

end #module
