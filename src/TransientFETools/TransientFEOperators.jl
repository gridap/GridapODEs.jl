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

function TransientFEOperator(trial::Union{FESpace,TransientTrialFESpace},
  test::FESpace,terms)
  # @santiagobadia : The time step value does not provide much info here...
  # so I can just use the one at 0.0, but not sure... we want to keep the
  # same assembler all the time...
  assem = SparseMatrixAssembler(test,trial_0)
  TransientFEOperatorFromTerms(trial,test,assem,terms)
end

# function (tfes::TransientFEOperatorFromTerms)(t::Real)
#   FEOperator
#   # trial(t),test,assem,terms
#   # return a fe_operator
# end
# see comment above, I don't think it has much sense

get_test(op::TransientFEOperatorFromTerms,t) = op.test

get_trial(op::TransientFEOperatorFromTerms,t) = op.trial(t)

function allocate_residual(op::TransientFEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  _, cellids = collect_cell_residual(uh,v,op.terms)
  # @santiagobadia : How do we do this? This assembler should be in the
  #
  allocate_vector(op.assem, cellids)
end

function residual!(b::AbstractVector,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  v = get_cell_basis(op.test)
  cellvecs, cellids = collect_cell_residual(t,uh,uh_t,v,op.terms)
  # @santiagobadia : not sure op will have Assembler
  assemble_vector!(b,op.assem, cellvecs, cellids)
  b
end

function allocate_jacobian(op::TransientFEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  # @santiagobadia : not sure op will have Assembler
  du = get_cell_basis(op.trial_0)
  # this is not a test function, it needs time... and it is not efficient
  v = get_cell_basis(op.test)
  _, cellidsrows, cellidscols = collect_cell_jacobian(uh,du,v,op.terms)
  allocate_matrix(op.assem, cellidsrows, cellidscols)
end

function jacobian!(A::AbstractMatrix,op::FEOperatorFromTerms,
  t::Real,uh,uh_t)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  du = get_cell_basis(op.trial_0)
  v = get_cell_basis(op.test)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian(t,uh,du,v,op.terms)
  assemble_matrix!(A,op.assem, cellmats, cellidsrows, cellidscols)
  A
end

function jacobian_t!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  du = get_cell_basis(op.trial_0)
  # This is not going to work, even though not needed Dir data here
  # @santiagobadia: get_trial can be expensive, how can we reuse its call?
  v = get_cell_basis(op.test)
  # to be implemented... collect_cell_jacobian_t
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian_t(t,uh,uh_t,du_t,v,op.terms)
  assemble_matrix!(A,op.assem, cellmats, cellidsrows, cellidscols)
  A
end

struct ODEOpFromFEOp <: ODEOperator
  feop::TransientFEOperator
end

function residual!(b::AbstractVector,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector)
  uhD = ...
  uhtD = ...
  uh = uhF, uhD
  uht = uhtF, uhtD
  # @santiagobadia : Here we need to reuse uh uht from somewhere, cache?
  # and create new FEFunctions by modifying the already assembled arrays
  residual!(b,op.feop,t,uh,uht)
end

# EvaluationFunction not FEFunction

function jacobian!(A::AbstractMatrix,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector)
  uhD = ...
  uhtD = ...
  uh = uhF, uhD
  uht = uhtF, uhtD
  # @santiagobadia : I want to do this only once!!!
  jacobian!(A,op.feop,t,uh,uht)
end

function jacobian_t!(A::AbstractMatrix,op::ODEOpFromFEOp,t,uh,uht)
  # @santiagobadia : Idem above
  jacobian_t!(A,op.feop,t,uh,uht)
end
