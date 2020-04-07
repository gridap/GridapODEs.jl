# const ∂t = time_derivative

"""
"""
abstract type TransientFEOperator <: GridapType end

(tfes::TransientFEOperator)(t::Real) = @notimplemented #::FEOperator

function get_test(op::TransientFEOperator)
  @abstractmethod
end

function get_trial(op::TransientFEOperator)
  @abstractmethod # time dependent
end

function allocate_residual(b::AbstractVector,op::TransientFEOperator,uh)
  @notimplemented
end

function allocate_jacobian(b::AbstractVector,op::TransientFEOperator,uh)
  @notimplemented
end

function residual!(b::AbstractVector,op::TransientFEOperator,t,uh,uht)
  @notimplemented
end

function jacobian!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht)
  @notimplemented
end

function jacobian_t!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht)
  @notimplemented
end

# @santiagobadia: missing allocates ...

struct TransientFEOperatorFromTerms <: TransientFEOperator
  trial::Union{FESpace,TransientTrialFESpace}
  trial_t::Union{FESpace,TransientTrialFESpace}
  test::FESpace
  assem_t::Function
  terms
  function TransientFEOperatorFromTerms(trial,test::FESpace,terms::FETerm...)
    new(trial,test,terms)
  end
end

function TransientFEOperator(trial::Union{FESpace,TransientTrialFESpace},
  test::FESpace,terms)
  TransientFEOperatorFromTerms(trial,∂t(trial),test,terms)
end

get_test(op::TransientFEOperatorFromTerms,t) = op.test

get_trial(op::TransientFEOperatorFromTerms,t) = op.trial(t)

function allocate_residual(op::TransientFEOperatorFromTerms,uh,assem)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  _, cellids = collect_cell_residual(uh,v,op.terms)
  allocate_vector(assem,cellids)
end

function residual!(b::AbstractVector,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,assem)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  v = get_cell_basis(op.test)
  cellvecs, cellids = collect_cell_residual(t,uh,uh_t,v,op.terms)
  assemble_vector!(b,assem,cellvecs,cellids)
  b
end

function allocate_jacobian(op::TransientFEOperatorFromTerms,uh,assem)
  @assert is_a_fe_function(uh)
  # @santiagobadia : not sure op will have Assembler
  du = get_cell_basis(op.trial_0)
  # this is not a test function, it needs time... and it is not efficient
  v = get_cell_basis(op.test)
  _, cellidsrows, cellidscols = collect_cell_jacobian(uh,du,v,op.terms)
  allocate_matrix(assem, cellidsrows, cellidscols)
end

function jacobian!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,assem)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  du = get_cell_basis(op.trial_0)
  v = get_cell_basis(op.test)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian(t,uh,du,v,op.terms)
  assemble_matrix!(A,assem, cellmats, cellidsrows, cellidscols)
  A
end

function jacobian_t!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,assem)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  du = get_cell_basis(op.trial_0)
  v = get_cell_basis(op.test)
  # to be implemented... collect_cell_jacobian_t
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian_t(t,uh,uh_t,du_t,v,op.terms)
  assemble_matrix!(A,assem, cellmats, cellidsrows, cellidscols)
  A
end

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
