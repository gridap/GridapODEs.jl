module TransientFETools

using Gridap.Helpers

export time_derivative
export ∂t

export TransientTrialFESpace
export get_trial_space
import Gridap.FESpaces: FESpace
import Gridap.FESpaces: SingleFieldFESpace
import Gridap.FESpaces: TrialFESpace
import Gridap.FESpaces: get_dirichlet_values

export TransientFETerm
import Gridap.FESpaces: FETerm
import Gridap.FESpaces: get_cell_jacobian
import Gridap.Geometry: Triangulation
import Gridap.Geometry: CellQuadrature

export TransientFEOperator
# export TransientFEOperatorFromTerms
# import Gridap.FESpaces: get_trial
using Gridap.FESpaces: Assembler
using Gridap.FESpaces: SparseMatrixAssembler
import GridapTimeStepper.ODETools: ODEOperator
import GridapTimeStepper.ODETools: residual!
import GridapTimeStepper.ODETools: jacobian!
import GridapTimeStepper.ODETools: jacobian_t!

export TransientFESolver
import Gridap.FESpaces: FESolver
import GridapTimeStepper.ODETools: ODESolver

include("TransientFESpaces.jl")

include("TransientFETerms.jl")

include("TransientFEOperators.jl")

include("TransientFESolvers.jl")

end #module
