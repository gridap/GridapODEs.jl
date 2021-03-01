module TransientFETools

using Test

using Gridap.Helpers

export ∂t

import GridapODEs.ODETools: ∂t, ∂tt
import GridapODEs.ODETools: time_derivative

export TransientTrialFESpace
export TransientMultiFieldFESpace
export test_transient_trial_fe_space
import Gridap.Fields: evaluate
import Gridap.Fields: evaluate!
import Gridap.MultiField: MultiFieldFESpace
using Gridap.FESpaces: FESpace
using Gridap.FESpaces: SingleFieldFESpace
using Gridap.FESpaces: TrialFESpace
using Gridap.FESpaces: ZeroMeanFESpace
using Gridap.FESpaces: get_free_dof_values
using Gridap.FESpaces: get_dirichlet_values
using Gridap.FESpaces: TrialFESpace!
using Gridap.FESpaces: HomogeneousTrialFESpace
using Gridap.FESpaces: jacobian

import Gridap.Geometry: Triangulation
import Gridap.CellData: Measure
using Gridap.FESpaces: ∫

export TransientFEOperator
export TransientAffineFEOperator
export TransientConstantFEOperator
using Gridap.FESpaces: Assembler
using Gridap.FESpaces: SparseMatrixAssembler
import GridapODEs.ODETools: allocate_cache
import GridapODEs.ODETools: update_cache!
import GridapODEs.ODETools: ODEOperator
import GridapODEs.ODETools: AffineODEOperator
import GridapODEs.ODETools: ConstantODEOperator
import GridapODEs.ODETools: allocate_residual
import GridapODEs.ODETools: allocate_jacobian
import GridapODEs.ODETools: residual!
import GridapODEs.ODETools: jacobian!
import GridapODEs.ODETools: jacobians!
import GridapODEs.ODETools: OperatorType
using GridapODEs.ODETools: Nonlinear
using GridapODEs.ODETools: Affine
using GridapODEs.ODETools: Constant
import Gridap.FESpaces: get_algebraic_operator
import Gridap.FESpaces: assemble_vector!
import Gridap.FESpaces: assemble_matrix_add!
import Gridap.FESpaces: allocate_vector
import Gridap.FESpaces: allocate_matrix
using Gridap.FESpaces: get_cell_shapefuns
using Gridap.FESpaces: get_cell_shapefuns_trial
using Gridap.FESpaces: collect_cell_vector
using Gridap.FESpaces: collect_cell_matrix
using Gridap.FESpaces: return_type
import Gridap.FESpaces: SparseMatrixAssembler
import Gridap.FESpaces: get_trial
import Gridap.FESpaces: get_test
using GridapODEs.ODETools: test_ode_operator
export test_transient_fe_operator

export TransientFESolver
import Gridap.FESpaces: FESolver
import GridapODEs.ODETools: ODESolver
import Gridap: solve
import Gridap.Algebra: solve!
import GridapODEs.ODETools: solve_step!
export test_transient_fe_solver

export TransientFEFunction
import Gridap.FESpaces: FEFunction
import Gridap.FESpaces: EvaluationFunction

export TransientFESolution
import Gridap: solve
import GridapODEs.ODETools: ODESolution
import GridapODEs.ODETools: GenericODESolution
import Base: iterate
export test_transient_fe_solution

export Transient2ndOrderFEOperator

include("TransientFESpaces.jl")

include("TransientFEOperators.jl")

include("ODEOperatorInterfaces.jl")

include("TransientFESolvers.jl")

include("TransientFESolutions.jl")

export FETerm
function FETerm(args...)
  Helpers.@unreachable """\n
  Function FETerm has been removed. The API for specifying the weak form has changed significantly.
  See the gridap/Tutorials repo for some examples of how to use the new API.
  This error message will be deleted in future versions.
  """
end

end #module
