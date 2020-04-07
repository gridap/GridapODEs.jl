# A concrete version of ODEOperator from a FEOperator
# struct ODEOperatorFromFEOperator <: ODEOperator
#   op::FEOperator
# end

struct ODEOpFromFEOp <: ODEOperator
  feop::TransientFEOperator
end

function allocate_state(op::ODEOpFromFEOp)
  Uh, Uht = _allocate_state(op,get_trial(op.feop))
  assem = op.assem_t(0.0)
  Uh, Uht, assem
end

function _allocate_state(fesp::FESpace)
  Uh = fesp
  Uht = HomogeneousTrialFESpace(fesp)
  Uh, Uht
end

function _allocate_state(fesp::TransientTrialFESpace)
  Uh = HomogeneousTrialFESpace(fesp.space)
  Uht = HomogeneousTrialFESpace(fesp.space)
  Uh, Uht
end

function update_state!(state,op::ODEOpFromFEOp,t::Real)
  _update_state!(state,get_trial(op.feop))
end

_update_state!(state,::FESpace) = nothing


function _update_state!(state,tfesp::TransientTrialFESpace)
  Uh, Uht, assem = state
  Uhnew = TrialFESpace!(get_dirichlet_values(Uh),Uh,tfesp.dirichlet_t(t))
  fun = tfesp.dirichlet_t
  fun_t = ∂t(fun)
  Uhtnew = TrialFESpace!(get_dirichlet_values(Uht),Uht,fun_t(t))
  Uhnew, Uhtnew, assem
end

function residual!(b::AbstractVector,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,op_state)
  Uh, Uht = op_state
  uh = FEFunction(uhF,Uh)
  uht = FEFunction(uhtF,Uht)
  residual!(b,op.feop,t,uh,uht)
end

function jacobian!(A::AbstractMatrix,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,op_state)
  Uh, Uht = op_state
  uh = FEFunction(uhF,Uh)
  uht = FEFunction(uhtF,Uht)
  jacobian!(A,op.feop,t,uh,uht)
end

function jacobian_t!(A::AbstractMatrix,op::ODEOpFromFEOp,t,uh,uht,op_state)
  Uh, Uht = op_state
  uh = FEFunction(uhF,Uh)
  uht = FEFunction(uhtF,Uht)
  jacobian_t!(A,op.feop,t,uh,uht)
end
