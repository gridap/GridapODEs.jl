module TransientFETools

using Gridap.Helpers

export time_derivative
export ∂t

export TransientTrialFESpace
export get_trial_space
import Gridap.FESpaces: FESpace
import Gridap.FESpaces: SingleFieldFESpace
import Gridap.FESpaces: TrialFESpace
import Gridap.FESpaces: TrialFESpace!
import Gridap.FESpaces: get_dirichlet_values

export TransientFETerm
import Gridap.FESpaces: FETerm
import Gridap.FESpaces: get_cell_jacobian
import Gridap.Geometry: Triangulation
import Gridap.Geometry: CellQuadrature
using Gridap.FESpaces: restrict
using Gridap.FESpaces: integrate
import Gridap.Geometry: get_cell_id

export TransientFEOperator
# export TransientFEOperatorFromTerms
# import Gridap.FESpaces: get_trial
using Gridap.FESpaces: Assembler
using Gridap.FESpaces: SparseMatrixAssembler
import GridapTimeStepper.ODETools: allocate_state
import GridapTimeStepper.ODETools: update_state!
import GridapTimeStepper.ODETools: ODEOperator
import GridapTimeStepper.ODETools: allocate_residual
import GridapTimeStepper.ODETools: allocate_jacobian
import GridapTimeStepper.ODETools: residual!
import GridapTimeStepper.ODETools: jacobian!
import GridapTimeStepper.ODETools: jacobian_t!
import Gridap.FESpaces: get_algebraic_operator
import Gridap.FESpaces: collect_cell_residual
import Gridap.FESpaces: collect_cell_jacobian
import Gridap.FESpaces: assemble_vector!
import Gridap.FESpaces: assemble_matrix!
import Gridap.FESpaces: alocate_vector
import Gridap.FESpaces: allocate_matrix
using Gridap.FESpaces: is_a_fe_function
using Gridap.FESpaces: is_a_fe_cell_basis
using Gridap.FESpaces: get_cell_basis

export TransientFESolver
import Gridap.FESpaces: FESolver
import GridapTimeStepper.ODETools: ODESolver
import Gridap.Algebra: solve!
import GridapTimeStepper.ODETools: solve_step!

export TransientFEFunction
import Gridap.FESpaces: FEFunction

export TransientFESolution
import Gridap: solve

include("TransientFESpaces.jl")

include("TransientFETerms.jl")

include("TransientFEOperators.jl")

include("ODEOperatorInterfaces.jl")

include("TransientFESolvers.jl")

include("TransientFESolutions.jl")

end #module
