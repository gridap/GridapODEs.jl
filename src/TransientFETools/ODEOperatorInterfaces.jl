# A concrete version of ODEOperator from a FEOperator
# struct ODEOperatorFromFEOperator <: ODEOperator
#   op::FEOperator
# end

struct ODEOpFromFEOp <: ODEOperator
  feop::TransientFEOperator
end

function allocate_state(op::ODEOpFromFEOp)
  Uh, Uht = _allocate_state(get_trial(op.feop))
  # assem = op.assem_t(0.0)
  assem = get_assembler(op.feop)
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
  _update_state!(state,get_trial(op.feop),t)
end

_update_state!(state,::FESpace,t) = nothing

function _update_state!(state,tfesp::TransientTrialFESpace,t::Real)
  Uh, Uht = state
  TrialFESpace!(Uh,0.0)
  TrialFESpace!(Uh,tfesp.dirichlet_t(t))
  fun = tfesp.dirichlet_t
  fun_t = ∂t(fun)
  TrialFESpace!(Uht,fun_t(t))
  Uh, Uht
end

function allocate_residual(op::ODEOpFromFEOp,uhF::AbstractVector,op_state)
  Uh, Uht = op_state
  uh = FEFunction(Uh,uhF)
  allocate_residual(op.feop,uh)#,op_state)
end

function allocate_jacobian(op::ODEOpFromFEOp,uhF::AbstractVector,op_state)
  Uh, Uht = op_state
  uh = FEFunction(Uh,uhF)
  allocate_jacobian(op.feop,uh,op_state)
end

function residual!(b::AbstractVector,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,op_state)
  Uh, Uht = op_state
  uh = FEFunction(Uh,uhF)
  uht = FEFunction(Uht,uhtF)
  residual!(b,op.feop,t,uh,uht,op_state)
end

function jacobian!(A::AbstractMatrix,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,op_state)
  Uh, Uht = op_state
  uh = FEFunction(Uh,uhF)
  uht = FEFunction(Uht,uhtF)
  jacobian!(A,op.feop,t,uh,uht,op_state)
end

function jacobian_t!(J::AbstractMatrix,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,dut_u::Real,op_state)
  Uh, Uht = op_state
  uh = FEFunction(Uh,uhF)
  uht = FEFunction(Uht,uhtF)
  jacobian_t!(J,op.feop,t,uh,uht,dut_u,op_state)
end
