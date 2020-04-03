const ∂t = time_derivative

# Do we need a TransientFEOperatorFromTerms or the only thing we need to do is
# to extend its interface. I would consider the second option
"""
"""
abstract type TransientFEOperator #<: FEOperator end
# @santiagobadia : ,: GridapType or FEOperator, it needs time to be a FEOperator
# thus probably not right subtyping...

function get_trial(op::TransientFEOperator,t)
  @abstractmethod
end

function residual!(b::AbstractVector,op::TransientFEOperator,uh,uht)
  @notimplemented
end

function jacobian!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht)
  @notimplemented
end

function jacobian_t!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht)
  @notimplemented
end

struct TransientFEOperatorFromTerms <: TransientFEOperator
  trial::FEspace
  test::FESpace
  assem::Assembler
  terms
  function TransientFEOperatorFromTerms(trial::FESpace,test::FESpace,assem::Assembler,terms::FETerm...)
    new(trial,test,assem,terms)
  end
end

function get_fe_operator(tfeop::TransientFEOperatorFromTerms,t::Real)
  # return a fe_operator
end

# TransientFEOperatorFromTerms(U::FESpace,V::FESpace,terms::TransientFETerm...)

# function TransientFEOperatorFromTerms(U::FESpace,U_t::FESpace,V::FESpace,terms::TransientFETerm...)
  # TransientFEOperatorFromTerms(V,t->U,t->U_t,terms)
# end

function get_trial(op::TransientFEOperatorFromTerms,t)
  get_trial(op.trial,t)
end


function residual!(b::AbstractVector,op::TransientFEOperatorFromTerms,uh,uht)
  # @santiagobadia : The uh and uht here are FEFunctions, created in the op
  # below. However, it is unclear to me whether we can keep the residual!
  # signature of FEOperator, since we need t to evaluate a transient operator
  # and both uh and uht... I have changed the method interface !!! But now it
  # should NOT subtype FEOperator !!!
  @notimplemented #yet
  A
end

function jacobian!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,t,uh,uht)
  @notimplemented #yet
  A
end

function jacobian_t!(A::AbstractMatrix,op::FEOperatorFromTerms,t,uh,uht)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  # This is not going to work, even though not needed Dir data here
  # @santiagobadia: get_trial can be expensive, how can we reuse its call?
  v = get_cell_basis(op.test)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian_unk_t(uh,uh_t,du_t,v,op.terms)
  assemble_matrix!(A,op.assem, cellmats, cellidsrows, cellidscols)
  A
end

function get_ode_operator(feop::TransientFEOperator)
  ODEOpFromFEOp(feop)
end

struct ODEOpFromFEOp <: ODEOperator
  feop::FEOperator
end

function residual!(b::AbstractVector,op::ODEOpFromFEOp,t,uh,uht)
  # @santiagobadia : uh are just free dof array, idem uht
  # Here we should think about strong bcs
  residual!(b,op.feop,t,uh,uht)
end

function jacobian!(A::AbstractMatrix,op::ODEOpFromFEOp,t,uh,uht)
  jacobian!(A,op.feop,t,uh,uht)
end

function jacobian_t!(A::AbstractMatrix,op::FEOperatorFromTerms,t,uh,uht)
  jacobian_t!(A,op.feop,t,uh,uht)
end
