module TransientFETools

using Test

using Gridap.Helpers

export time_derivative
export ∂t

export TransientTrialFESpace
export test_transient_trial_fe_space
import Gridap.Fields: evaluate
import Gridap.Fields: evaluate!
import Gridap.MultiField: MultiFieldFESpace
using Gridap.FESpaces: FESpace
using Gridap.FESpaces: SingleFieldFESpace
using Gridap.FESpaces: TrialFESpace
using Gridap.FESpaces: get_free_values
using Gridap.FESpaces: get_dirichlet_values
using Gridap.FESpaces: TrialFESpace!
using Gridap.FESpaces: HomogeneousTrialFESpace

export TransientFETerm
import Gridap.FESpaces: FETerm
import Gridap.FESpaces: get_cell_jacobian
import Gridap.Geometry: Triangulation
import Gridap.Geometry: CellQuadrature
using Gridap.FESpaces: restrict
using Gridap.FESpaces: integrate
import Gridap.Geometry: get_cell_id

export TransientFEOperator
using Gridap.FESpaces: Assembler
using Gridap.FESpaces: SparseMatrixAssembler
import GridapTimeStepper.ODETools: allocate_cache
import GridapTimeStepper.ODETools: update_cache!
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
import Gridap.FESpaces: assemble_matrix_add!
import Gridap.FESpaces: allocate_vector
import Gridap.FESpaces: allocate_matrix
using Gridap.FESpaces: is_a_fe_function
using Gridap.FESpaces: is_a_fe_cell_basis
using Gridap.FESpaces: get_cell_basis
using GridapTimeStepper.ODETools: test_ode_operator
export test_transient_fe_operator

export TransientFESolver
import Gridap.FESpaces: FESolver
import GridapTimeStepper.ODETools: ODESolver
import Gridap: solve
import Gridap.Algebra: solve!
import GridapTimeStepper.ODETools: solve_step!
export test_transient_fe_solver

export TransientFEFunction
import Gridap.FESpaces: FEFunction

export TransientFESolution
import Gridap: solve
import GridapTimeStepper.ODETools: ODESolution
import GridapTimeStepper.ODETools: GenericODESolution
import Base: iterate
export test_transient_fe_solution

include("TransientFESpaces.jl")

include("TransientFETerms.jl")

include("TransientFEOperators.jl")

include("ODEOperatorInterfaces.jl")

include("TransientFESolvers.jl")

include("TransientFESolutions.jl")

end #module
