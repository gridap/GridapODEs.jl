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
using Gridap.FESpaces: get_dirichlet_dof_values
using Gridap.FESpaces: TrialFESpace!
using Gridap.FESpaces: HomogeneousTrialFESpace
using Gridap.FESpaces: jacobian

import Gridap.Geometry: Triangulation
import Gridap.CellData: Measure
using Gridap.FESpaces: ∫

export TransientFEOperator
export TransientAffineFEOperator
export TransientConstantFEOperator
export TransientConstantMatrixFEOperator
using Gridap.FESpaces: Assembler
using Gridap.FESpaces: SparseMatrixAssembler
import GridapODEs.ODETools: allocate_cache
import GridapODEs.ODETools: update_cache!
import GridapODEs.ODETools: ODEOperator
import GridapODEs.ODETools: AffineODEOperator
import GridapODEs.ODETools: ConstantODEOperator
import GridapODEs.ODETools: ConstantMatrixODEOperator
import GridapODEs.ODETools: allocate_residual
import GridapODEs.ODETools: allocate_jacobian
import GridapODEs.ODETools: residual!
import GridapODEs.ODETools: jacobian!
import GridapODEs.ODETools: jacobians!
import GridapODEs.ODETools: OperatorType
using GridapODEs.ODETools: Nonlinear
using GridapODEs.ODETools: Affine
using GridapODEs.ODETools: Constant
using GridapODEs.ODETools: ConstantMatrix
import Gridap.FESpaces: get_algebraic_operator
import Gridap.FESpaces: assemble_vector!
import Gridap.FESpaces: assemble_matrix_add!
import Gridap.FESpaces: allocate_vector
import Gridap.FESpaces: allocate_matrix
using Gridap.FESpaces: get_fe_basis
using Gridap.FESpaces: get_trial_fe_basis
using Gridap.FESpaces: collect_cell_vector
using Gridap.FESpaces: collect_cell_matrix
using Gridap.FESpaces: return_type
import Gridap.FESpaces: SparseMatrixAssembler
import Gridap.FESpaces: get_trial
import Gridap.FESpaces: get_test
using GridapODEs.ODETools: test_ode_operator
export test_transient_fe_operator

import Gridap.FESpaces: FESolver
import GridapODEs.ODETools: ODESolver
import Gridap: solve
import Gridap.Algebra: solve!
import GridapODEs.ODETools: solve_step!
export test_transient_fe_solver

export TransientFEFunction
import Gridap.FESpaces: FEFunction
import Gridap.FESpaces: SingleFieldFEFunction
import Gridap.FESpaces: EvaluationFunction
import Gridap.MultiField: MultiFieldFEFunction
import Gridap.MultiField: num_fields

export TransientFESolution
import Gridap: solve
import GridapODEs.ODETools: ODESolution
import GridapODEs.ODETools: GenericODESolution
import Base: iterate
export test_transient_fe_solution

export TransientCellField
using Gridap.CellData: CellField
using Gridap.CellData: CellFieldAt
using Gridap.CellData: GenericCellField
using Gridap.MultiField: MultiFieldCellField
using Gridap.FESpaces: FEBasis
import Gridap.CellData: get_data
import Gridap.CellData: get_triangulation
import Gridap.CellData: DomainStyle
import Gridap.CellData: gradient
import Gridap.CellData: ∇∇
import Gridap.CellData: change_domain
import Gridap.FESpaces: BasisStyle

using PartitionedArrays: get_part, map_parts, get_part_ids
using PartitionedArrays: PSparseMatrix
using GridapDistributed: local_views
using GridapDistributed: DistributedCellDatum
using GridapDistributed: DistributedMeasure
using GridapDistributed: DistributedCellField
using GridapDistributed: DistributedSingleFieldFESpace
using GridapDistributed: DistributedMultiFieldFESpace
using GridapDistributed: DistributedSingleFieldFEFunction
using GridapDistributed: DistributedMultiFieldFEFunction
using Gridap.Helpers
using Gridap.Fields
using Gridap.Arrays
using Gridap.CellData
using Gridap.Geometry: SkeletonPair

import Gridap.TensorValues: inner, outer, double_contraction, symmetric_part
import LinearAlgebra: det, tr, cross, dot, ⋅
import Base: inv, abs, abs2, *, +, -, /, adjoint, transpose, real, imag, conj

export TransientDistributedCellField
export TransientSingleFieldDistributedCellField

include("TransientFESpaces.jl")

include("TransientCellField.jl")

include("TransientMultiFieldCellField.jl")

include("TransientDistributedCellField.jl")

include("TransientFEOperators.jl")

include("ODEOperatorInterfaces.jl")

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
