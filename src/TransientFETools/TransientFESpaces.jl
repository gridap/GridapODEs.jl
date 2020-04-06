struct TransientTrialFESpace
  space::FESpace
  dirichlet_t::Vector{<:Function}
  dirichlet_values::AbstractVector
end

function TransientTrialFESpace(U::FESpace,dir_fun)
 dirichlet_vals = similar(get_dirichlet_values(U))
  TransientTrialFESpace(U,dir_fun,dirichlet_values)
end

function (fes::TransientTrialFESpace)(t::Real)
  TrialFESpace!(fes.dirichlet_values,fes.space,fes.dirichlet_t(t))
end

(tfes::FESpace)(t::Real) = tfes

∂t(fes::TransientTrialFESpace) = TransientTrialFESpace(fes.space,∂t.(fes.dirichlet_t))

∂t(fes::TrialFESpace) = HomogeneousTrialFESpace(fes.space)


# to gridap
function HomogeneousTrialFESpace(U::FESpace)
 dirichlet_vals = similar(get_dirichlet_values(U))
  fill!(dirichlet_vals,zero(eltype(dirichlet_vals)))
  # to be done in Gridap
  # TrialFESpace with dirichletvalues
  TrialFESpace(U,dirichlet_values)
end
