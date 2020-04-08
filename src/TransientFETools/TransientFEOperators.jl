# const ∂t = time_derivative

"""
"""
abstract type TransientFEOperator <: GridapType end

# @santiagobadia : To be eliminated ?
(tfes::TransientFEOperator)(t::Real) = @notimplemented #::FEOperator

function get_test(op::TransientFEOperator)
  @abstractmethod
end

function get_trial(op::TransientFEOperator)
  @abstractmethod # time dependent
end

function allocate_residual(op::TransientFEOperator,uh)#,state)
  @notimplemented
end

function allocate_jacobian(op::TransientFEOperator,uh,state)
  @notimplemented
end

function residual!(b::AbstractVector,op::TransientFEOperator,t,uh,uht,state)
  @notimplemented
end

function jacobian!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht,state)
  @notimplemented
end

function jacobian_t!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht,duht_du,state)
  @notimplemented
end

function get_algebraic_operator(feop::TransientFEOperator)
  ODEOpFromFEOp(feop)
end

get_assembler(feop::TransientFEOperator) = @notimplemented

# @santiagobadia: missing allocates ...

struct TransientFEOperatorFromTerms <: TransientFEOperator
  trial::Union{FESpace,TransientTrialFESpace}
  trial_t::Union{FESpace,TransientTrialFESpace}
  test::FESpace
  assem_t::Assembler
  terms
end

get_assembler(feop::TransientFEOperatorFromTerms) = feop.assem_t

function TransientFEOperator(trial::Union{FESpace,TransientTrialFESpace},
  test::FESpace,terms::Union{FETerm,TransientFETerm}...)
  # @santiagobadia : I am here assem_t ... can we create assem here?
  assem_t = SparseMatrixAssembler(test,trial(0.0))
  TransientFEOperatorFromTerms(trial,∂t(trial),test,assem_t,terms...)
end

get_test(op::TransientFEOperatorFromTerms) = op.test

get_trial(op::TransientFEOperatorFromTerms) = op.trial

get_trial(op::TransientFEOperatorFromTerms,t) = op.trial(t)

function allocate_residual(op::TransientFEOperatorFromTerms,uh)#,state)#,assem)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  _, cellids =
  collect_cell_residual(0.0,uh,uh,v,op.terms)
  allocate_vector(op.assem_t,cellids)
end

function residual!(b::AbstractVector,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,state)#,assem)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  v = get_cell_basis(op.test)
  cellvecs, cellids = collect_cell_residual(t,uh,uh_t,v,op.terms)
  assemble_vector!(b,op.assem_t,cellvecs,cellids)
  b
end

function allocate_jacobian(op::TransientFEOperatorFromTerms,uh,state)
  Uh, Uht = state
  @assert is_a_fe_function(uh)
  # @santiagobadia : not sure op will have Assembler
  du = get_cell_basis(Uh)
  # this is not a test function, it needs time... and it is not efficient
  v = get_cell_basis(op.test)
  _, cellidsrows, cellidscols = collect_cell_jacobian(0.0,uh,uh,du,v,op.terms)
  allocate_matrix(op.assem_t, cellidsrows, cellidscols)
end

function jacobian!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,state)#,assem)
  Uh, Uht = state
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  du = get_cell_basis(Uh)
  v = get_cell_basis(op.test)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian(t,uh,uh_t,du,v,op.terms)
  assemble_matrix!(A,op.assem_t, cellmats, cellidsrows, cellidscols)
  A
end

function jacobian_t!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,duht_du::Real,state)#,assem)
  Uh, Uht = state
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  du_t = get_cell_basis(Uh)
  v = get_cell_basis(op.test)
  # to be implemented... collect_cell_jacobian_t
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian_t(t,uh,uh_t,du_t,v,op.terms)
  assemble_matrix!(A,op.assem_t, cellmats, cellidsrows, cellidscols)
  A
end

function test_transient_fe_operator(op::TransientFEOperator,uh)
  odeop = get_algebraic_operator(op)
  @test isa(odeop,ODEOperator)
  state = allocate_state(odeop)
  V = get_test(op)
  @test isa(V,FESpace)
  U = get_trial(op)
  @test isa(U,Union{FESpace,TransientTrialFESpace})
  U0 = U(0.0)
  @test isa(U0,FESpace)
  # uh = FEFunction(U0,0.0)
  r = allocate_residual(op,uh)
  @test isa(r,AbstractVector)
  residual!(r,op,0.0,uh,uh,state)
  # r2 = residual(op)
  @test isa(r,AbstractVector)
  J = allocate_jacobian(op,uh,state)
  @test isa(J,AbstractMatrix)
  jacobian!(J,op,0.0,uh,uh,state)
  @test isa(J,AbstractMatrix)
  jacobian_t!(J,op,0.0,uh,uh,1.0,state)
  @test isa(J,AbstractMatrix)
  # uhF = get_free_values(uh)
  # test_ode_operator(odeop,0.0,uhF,uhF)
  # get_assembler(feop::TransientFEOperator) = @notimplemented
  true
end
