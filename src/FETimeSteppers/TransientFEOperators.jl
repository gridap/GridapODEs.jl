# Do we need a TransientFEOperatorFromTerms or the only thing we need to do is
# to extend its interface. I would consider the second option
"""
"""
abstract type TransientFEOperator <: GridapType end

struct TransientFEOperatorFromTerms <: TransientFEOperator
  test::FESpace
  trial::Function
  trial_t::Function
  terms::Tuple
end

TransientFEOperatorFromTerms(U::FESpace,V::FESpace,terms::TransientFETerm...)

function TransientFEOperatorFromTerms(U::FESpace,U_t::FESpace,V::FESpace,terms::TransientFETerm...)
  TransientFEOperatorFromTerms(V,t->U,t->U_t,terms)
end

function jacobian_unk_t!(A::AbstractMatrix,op::FEOperator,uh)
  @notimplemented
  A
end

function jacobian_unk_t!(A::AbstractMatrix,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian_unk_t(uh,uh_t,du_t,v,op.terms)
  assemble_matrix!(A,op.assem, cellmats, cellidsrows, cellidscols)
  A
end

function get_ode_operator(feop::ODEOperator)
  ODEOpFromFEOp(feop)
end

struct ODEOpFromFEOp <: ODEOperator
  feop::FEOperator
end

