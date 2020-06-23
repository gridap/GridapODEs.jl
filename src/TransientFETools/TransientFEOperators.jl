"""
A transient version of the `Gridap` `FEOperator` that depends on time
"""
abstract type TransientFEOperator <: GridapType end

"""
Returns the test space
"""
function get_test(op::TransientFEOperator)
  @abstractmethod
end

"""
Returns the (possibly) time-dependent trial space
"""
function get_trial(op::TransientFEOperator)
  @abstractmethod # time dependent
end

function allocate_residual(op::TransientFEOperator,uh,cache)
  @abstractmethod
end

function allocate_jacobian(op::TransientFEOperator,uh,cache)
  @notimplemented
end

"""
Idem as `residual!` of `ODEOperator`
"""
function residual!(b::AbstractVector,op::TransientFEOperator,t,uh,uht,cache)
  @abstractmethod
end

"""
Idem as `residual!` of `ODEOperator`
"""
function jacobian!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht,cache)
  @abstractmethod
end

"""
Idem as `jacobian_t!` of `ODEOperator`
"""
function jacobian_t!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht,duht_du,cache)
  @abstractmethod
end

"""
Returns the assembler, which is constant for all time steps for a given FE
operator.

Note: adaptive FE spaces involve to generate new FE spaces and
corresponding operators, due to the ummutable approach in `Gridap`
"""
get_assembler(feop::TransientFEOperator) = @abstractmethod


# Default API

"""
Returns a `ODEOperator` wrapper of the `TransientFEOperator` that can be
straightforwardly used with the `ODETools` module.
"""
get_algebraic_operator(feop::TransientFEOperator) = @abstractmethod

# @fverdugo This function is just in case we need to override it in the future for some specialization.
# This default implementation is enough for the moment.
function allocate_cache(op::TransientFEOperator)
  nothing
end

function update_cache!(cache::Nothing,op::TransientFEOperator,t::Real)
  nothing
end

# Specializations

"""
Transient FE operator that is defined by a set of `TransientFETerm` (or `FETerm`)
"""
struct TransientFEOperatorFromTerms <: TransientFEOperator
  trial
  trial_t
  test::FESpace
  assem_t::Assembler
  terms
  type::OperatorType
  function TransientFEOperatorFromTerms(trial,trial_t,test::FESpace,assem::Assembler,terms...)
    new(trial,trial_t,test,assem,Nonlinear(),terms...)
  end
end

function TransientConstantFEOperator(trial,test,terms...)
  assem_t = SparseMatrixAssembler(test,evaluate(trial,nothing))
  TransientFEOperatorFromTerms(trial,∂t(trial),test,assem_t,Constant(),terms...)
end

function TransientAffineFEOperator(trial,test,terms...)
  assem_t = SparseMatrixAssembler(test,evaluate(trial,nothing))
  TransientFEOperatorFromTerms(trial,∂t(trial),test,assem_t,Affine(),terms...)
end

function TransientFEOperator(trial,test,terms...)
  assem_t = SparseMatrixAssembler(test,evaluate(trial,nothing))
  # if (type == "nonlinear")
  #   optype = Nonlinear()
  # elseif (type == "affine")
  #   optype = Affine()
  # elseif (type == "constant")
  #   optype = Constant()
  # else
  #   error("Operator type not defined, it must be nonlinear, affine or constant")
  # end
  TransientFEOperatorFromTerms(trial,∂t(trial),test,assem_t,Nonlinear(),terms...)
end

get_assembler(feop::TransientFEOperatorFromTerms) = feop.assem_t

get_test(op::TransientFEOperatorFromTerms) = op.test

get_trial(op::TransientFEOperatorFromTerms) = op.trial

function allocate_residual(op::TransientFEOperatorFromTerms,uh,cache)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  vecdata = collect_cell_residual(0.0,uh,uh,v,op.terms)
  allocate_vector(op.assem_t,vecdata)
end

function residual!(b::AbstractVector,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,cache)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  v = get_cell_basis(op.test)
  vecdata = collect_cell_residual(t,uh,uh_t,v,op.terms)
  assemble_vector!(b,op.assem_t,vecdata)
  b
end

function allocate_jacobian(op::TransientFEOperatorFromTerms,uh,cache)
  Uh = evaluate(op.trial,nothing)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(Uh)
  v = get_cell_basis(op.test)
  matdata = collect_cell_jacobian(0.0,uh,uh,du,v,op.terms)
  allocate_matrix(op.assem_t,matdata)
end

function jacobian!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,cache)
  Uh = evaluate(op.trial,nothing)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  du = get_cell_basis(Uh)
  v = get_cell_basis(op.test)
  matdata = collect_cell_jacobian(t,uh,uh_t,du,v,op.terms)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function jacobian_t!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,duht_du::Real,cache)
  Uh = evaluate(op.trial,nothing)
  @assert is_a_fe_function(uh_t)
  du_t = get_cell_basis(Uh)
  v = get_cell_basis(op.test)
  matdata = collect_cell_jacobian_t(t,uh,uh_t,du_t,v,duht_du,op.terms)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function get_algebraic_operator(feop::TransientFEOperatorFromTerms)
  ODEOpFromFEOp{typeof(feop.type)}(feop)
end

# Tester

function test_transient_fe_operator(op::TransientFEOperator,uh)
  odeop = get_algebraic_operator(op)
  @test isa(odeop,ODEOperator)
  cache = allocate_cache(op)
  V = get_test(op)
  @test isa(V,FESpace)
  U = get_trial(op)
  U0 = U(0.0)
  @test isa(U0,FESpace)
  r = allocate_residual(op,uh,cache)
  @test isa(r,AbstractVector)
  residual!(r,op,0.0,uh,uh,cache)
  @test isa(r,AbstractVector)
  J = allocate_jacobian(op,uh,cache)
  @test isa(J,AbstractMatrix)
  jacobian!(J,op,0.0,uh,uh,cache)
  @test isa(J,AbstractMatrix)
  jacobian_t!(J,op,0.0,uh,uh,1.0,cache)
  @test isa(J,AbstractMatrix)
  cache = update_cache!(cache,op,0.0)
  true
end
