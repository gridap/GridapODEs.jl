struct TransientTrialFESpace # not sure subtyping...
  space::SingleFieldFESpace
  dirichlet_t::Function
end

function get_trial(tfes::TransientTrialFESpace,t::Real)
  TrialFESpace(tfes.space,dirichlet_t(t))
end

# and here we extend the FESpace interface
get_trial(fesp::TrialFESpace,t::Real) = fesp
