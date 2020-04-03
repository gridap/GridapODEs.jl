struct TransientTrialFESpace
  space::SingleFieldFESpace
  dirichlet_t::Function
end

(tfes::TransientTrialFESpace)(t::Real) = TrialFESpace(tfes.space,dirichlet_t(t))

(tfes::FESpace)(t::Real) = tfes

∂t(fes::TransientTrialFESpace) = TransientTrialFESpace(fes.space,∂t(fes.dirichlet_t))

∂t(fes::TrialFESpace) = HomogeneousTrialFESpace(fes.space)

# to gridap
function HomogeneousTrialFESpace(U::FESpace)
 dirichlet_vals = similar(get_dirichlet_values(U))
  fill!(dirichlet_vals,zero(eltype(dirichlet_vals)))
  # to be done in Gridap
  TrialFESpace(U,dirichlet_values)
end
