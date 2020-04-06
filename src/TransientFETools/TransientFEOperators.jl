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
  trial_0::FESpace
  test::FESpace
  # @santiagobadia : not sure it should be here
  assem::Assembler
  terms
  function TransientFEOperatorFromTerms(trial,test::FESpace,assem::Assembler,terms::FETerm...)
    new(trial,test,assem,terms)
  end
end

function TransientFEOperator(trial::Union{FESpace,TransientTrialFESpace},
  test::FESpace,terms)
  # @santiagobadia : The time step value does not provide much info here...
  # so I can just use the one at 0.0, but not sure... we want to keep the
  # same assembler all the time...
  tr0 = trial(0.0)
  assem = SparseMatrixAssembler(test,tr0)
  TransientFEOperatorFromTerms(trial,tr0,test,assem,terms)
end

function (tfes::TransientFEOperatorFromTerms)(t::Real)
  FEOperator
  # trial(t),test,assem,terms
  # return a fe_operator
end
# see comment above, I don't think it has much sense

get_test(op::TransientFEOperatorFromTerms,t) = op.test

get_trial(op::TransientFEOperatorFromTerms,t) = op.trial(t)

function residual!(b::AbstractVector,op::TransientFEOperatorFromTerms,t,uh,uht)
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

function jacobian_t!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,t,uh,uht)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  # This is not going to work, even though not needed Dir data here
  # @santiagobadia: get_trial can be expensive, how can we reuse its call?
  v = get_cell_basis(op.test)
  # to be implemented... collect_cell_jacobian_t
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian_t(uh,uh_t,du_t,v,op.terms)
  assemble_matrix!(A,op.assem, cellmats, cellidsrows, cellidscols)
  A
end

function get_ode_operator(feop::TransientFEOperator)
  ODEOpFromFEOp(feop)
end

struct ODEOpFromFEOp <: ODEOperator
  feop::TransientFEOperator
end

function residual!(b::AbstractVector,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector)
  # @santiagobadia : uh are just free dof array, idem uht
  # Here we should think about strong bcs
  residual!(b,op.feop,t,uh,uht)
end

# EvaluationFunction not FEFunction

function jacobian!(A::AbstractMatrix,op::ODEOpFromFEOp,t,uh,uht)
  jacobian!(A,op.feop,t,uh,uht)
end

function jacobian_t!(A::AbstractMatrix,op::ODEOpFromFEOp,t,uh,uht)
  jacobian_t!(A,op.feop,t,uh,uht)
end
