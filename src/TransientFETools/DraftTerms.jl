
struct TransientFEOperatorFromTerms <: TransientFEOperator
  trial::Union{FESpace,TransientTrialFESpace}
  test::FESpace
  # assem::Assembler
  terms
  function TransientFEOperatorFromTerms(trial,test::FESpace,terms::FETerm...)
    new(trial,test,terms)
    # new(trial,test,assem,terms)
  end
end

function TransientFEOperator(trial::Union{FESpace,TransientTrialFESpace},
  test::FESpace,terms)
  # @santiagobadia : The time step value does not provide much info here...
  # so I can just use the one at 0.0, but not sure... we want to keep the
  # same assembler all the time...
  assem = SparseMatrixAssembler(test,trial(0.0))
  TransientFEOperatorFromTerms(trial,test,assem,terms)
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

# """
# """
# function FEOperator(trial::FESpace,test::FESpace,terms::FETerm...)
#   assem = SparseMatrixAssembler(test,trial)
#   FEOperator(trial,test,assem,terms...)
# end

function get_test(op::FEOperatorFromTerms)
  op.test
end

function get_trial(op::FEOperatorFromTerms)
  op.trial
end

function allocate_residual(op::TransientFEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  _, cellids = collect_cell_residual(uh,v,op.terms)
  allocate_vector(op.assem, cellids)
end

function residual!(b::AbstractVector,op::TransientFEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  cellvecs, cellids = collect_cell_residual(uh,v,op.terms)
  assemble_vector!(b,op.assem, cellvecs, cellids)
  b
end

function allocate_jacobian(op::TransientFEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  _, cellidsrows, cellidscols = collect_cell_jacobian(uh,du,v,op.terms)
  allocate_matrix(op.assem, cellidsrows, cellidscols)
end

function jacobian!(A::AbstractMatrix,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian(uh,du,v,op.terms)
  assemble_matrix!(A,op.assem, cellmats, cellidsrows, cellidscols)
  A
end

function residual_and_jacobian!(b::AbstractVector,A::AbstractMatrix,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  data = collect_cell_jacobian_and_residual(uh,du,v,op.terms)
  assemble_matrix_and_vector!(A, b, op.assem,data...)
  (b,A)
end

function residual_and_jacobian(op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  data = collect_cell_jacobian_and_residual(uh,du,v,op.terms)
  A, b = assemble_matrix_and_vector(op.assem,data...)
  (b, A)
end
