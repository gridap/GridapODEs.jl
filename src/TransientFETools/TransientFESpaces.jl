struct TransientTrialFESpace
  space::FESpace
  dirichlet_t#::Vector{<:Function}
  dirichlet_values::AbstractVector
end

# function TransientTrialFESpace(U::FESpace,dir_fun)
  # TransientTrialFESpace(U,[dir_fun])
# end

function TransientTrialFESpace(U::FESpace,objects)
 dirichlet_vals = similar(get_dirichlet_values(U))
  TransientTrialFESpace(U,objects,dirichlet_vals)
end

function (fes::TransientTrialFESpace)(t::Real)
  # dir_t = [ud(t) for ud in fes.dirichlet_t]
  dir_t = fes.dirichlet_t(t)
  TrialFESpace(fes.dirichlet_values,fes.space,dir_t)
end

# # @santiagobadia : Using private methods from Gridal
# using Gridap.FESpaces: _prepare_trial_cell_basis
# using Gridap.FESpaces: constraint_style
#
# function TrialFESpace(dir_values::AbstractVector,fesp,dir_funs)
#   TrialFESpace(fesp,dir_values,_prepare_trial_cell_basis(fesp),constraint_style(fesp))
# end


(tfes::FESpace)(t::Real) = tfes

∂t(fes::TransientTrialFESpace) = TransientTrialFESpace(fes.space,∂t.(fes.dirichlet_t))

∂t(fes::TrialFESpace) = HomogeneousTrialFESpace(fes.space)


# to gridap
function HomogeneousTrialFESpace(U::FESpace)
  # @santiagobadia : To be improved
  TrialFESpace(U,0.0)
 # dirichlet_vals = similar(get_dirichlet_values(U))
  # fill!(dirichlet_vals,zero(eltype(dirichlet_vals)))
  # to be done in Gridap
  # TrialFESpace with dirichletvalues
  # TrialFESpace(U,dirichlet_vals)
end

function test_transient_trial_fe_space(Uh::TransientTrialFESpace)
  Uh0=Uh(0.0)
  Uh00=Uh0(0.0)
  Uht=∂t(Uh)
  Uht0=Uht(0.0)
  Uh0t = ∂t(Uh0)
  true
end
