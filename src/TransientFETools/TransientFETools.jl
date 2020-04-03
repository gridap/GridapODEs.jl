module TransientFETools

using Gridap.Helpers

export time_derivative
export ∂t

export TransientTrialFESpace
export get_trial_space
import Gridap.FESpaces: FESpace
import Gridap.FESpaces: SingleFieldFESpace
import Gridap.FESpaces: TrialFESpace

export TransientFETerm
import Gridap.FESpaces: FETerm
import Gridap.FESpaces: get_cell_jacobian
import Gridap.Geometry: Triangulation
import Gridap.Geometry: CellQuadrature

export TransientFEOperator
export TransientFEOperatorFromTerms
export TransientFEOperatorFromTerms
import Gridap.FESpaces: get_trial
import GridapTimeStepper.ODETools: ODEOperator
import GridapTimeStepper.ODETools: residual!
import GridapTimeStepper.ODETools: jacobian!
import GridapTimeStepper.ODETools: jacobian_t!

include("TransientFESpaces.jl")

include("TransientFETerms.jl")

end #module
