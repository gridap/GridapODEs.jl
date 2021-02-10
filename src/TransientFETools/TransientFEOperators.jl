"""
A transient version of the `Gridap` `FEOperator` that depends on time
"""
abstract type TransientFEOperator{C<:OperatorType} <: GridapType end

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
Idem as `jacobian_and_jacobian_t!` of `ODEOperator`
"""
function jacobian_and_jacobian_t!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht,duht_du,cache)
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
function get_algebraic_operator(feop::TransientFEOperator{C}) where C
  ODEOpFromFEOp{C}(feop)
end

OperatorType(::Type{<:TransientFEOperator{C}}) where C = C

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
Transient FE operator that is defined by a transient Weak form
"""
struct TransientFEOperatorFromWeakForm{C} <: TransientFEOperator{C}
  res::Function
  jac::Function
  jac_t::Function
  assem_t::Assembler
  trial
  trial_t
  test::FESpace
end

function TransientConstantFEOperator(m::Function,a::Function,b::Function,
  trial,test)
  res(t,u,ut,v) = m(ut,v) + a(u,v) - b(v)
  jac(t,u,ut,du,v) = a(du,v)
  jac_t(t,u,ut,dut,v) = m(dut,v)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{Constant}(res,jac,jac_t,assem_t,trial,∂t(trial),test)
end

function TransientAffineFEOperator(m::Function,a::Function,b::Function,
  trial,test)
  res(t,u,ut,v) = m(t,ut,v) + a(t,u,v) - b(t,v)
  jac(t,u,ut,du,v) = a(t,du,v)
  jac_t(t,u,ut,dut,v) = m(t,dut,v)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{Affine}(res,jac,jac_t,assem_t,trial,∂t(trial),test)
end

function TransientFEOperator(res::Function,jac::Function,jac_t::Function,
  trial,test)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{Nonlinear}(res,jac,jac_t,assem_t,trial,∂t(trial),test)
end
function TransientFEOperator(res::Function,trial,test)
  jac(t,u,ut,du,dv) = jacobian(x->res(t,x,ut,dv),u)
  jac_t(t,u,ut,dut,dv) = jacobian(xt->res(t,u,xt,dv),ut)
  TransientFEOperator(res,jac,jac_t,trial,test)
end

function SparseMatrixAssembler(
  trial::Union{TransientTrialFESpace,TransientMultiFieldTrialFESpace},
  test::FESpace)
  SparseMatrixAssembler(evaluate(trial,nothing),test)
end

get_assembler(op::TransientFEOperatorFromWeakForm) = op.assem_t

get_test(op::TransientFEOperatorFromWeakForm) = op.test

get_trial(op::TransientFEOperatorFromWeakForm) = op.trial

function allocate_residual(op::TransientFEOperatorFromWeakForm,uh::FEFunction,cache)
  v = get_cell_shapefuns(get_test(op))
  vecdata = collect_cell_vector(op.res(0.0,uh,uh,v))
  allocate_vector(op.assem_t,vecdata)
end

function residual!(b::AbstractVector,op::TransientFEOperatorFromWeakForm,
  t::Real,uh::FEFunction,uh_t::FEFunction,cache)
  v = get_cell_shapefuns(get_test(op))
  vecdata = collect_cell_vector(op.res(t,uh,uh_t,v))
  assemble_vector!(b,op.assem_t,vecdata)
  b
end

function allocate_jacobian(op::TransientFEOperatorFromWeakForm,uh::FEFunction,cache)
  matdata_j = matdata_jacobian(op,0.0,uh,uh)
  matdata_jt = matdata_jacobian_t(op,0.0,uh,uh,0.0)
  matdata = vcat_matdata(matdata_j,matdata_jt)
  allocate_matrix(op.assem_t,matdata)
end

function jacobian!(A::AbstractMatrix,op::TransientFEOperatorFromWeakForm,
  t::Real,uh::FEFunction,uh_t::FEFunction,cache)
  matdata = matdata_jacobian(op,t,uh,uh_t)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function jacobian_t!(A::AbstractMatrix,op::TransientFEOperatorFromWeakForm,
  t::Real,uh,uh_t,duht_du::Real,cache)
  matdata = matdata_jacobian_t(op,t,uh,uh_t,duht_du)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function jacobian_and_jacobian_t!(A::AbstractMatrix,op::TransientFEOperatorFromWeakForm,
  t::Real,uh,uh_t,duht_du::Real,cache)
  matdata_j = matdata_jacobian(op,t,uh,uh_t)
  matdata_jt = matdata_jacobian_t(op,t,uh,uh_t,duht_du)
  matdata = vcat_matdata(matdata_j,matdata_jt)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function vcat_matdata(matdata_j,matdata_jt)
  term_to_cellmat_j, term_to_cellidsrows_j, term_to_cellidscols_j = matdata_j
  term_to_cellmat_jt, term_to_cellidsrows_jt, term_to_cellidscols_jt = matdata_jt

  term_to_cellmat = vcat(term_to_cellmat_j,term_to_cellmat_jt)
  term_to_cellidsrows = vcat(term_to_cellidsrows_j,term_to_cellidsrows_jt)
  term_to_cellidscols = vcat(term_to_cellidscols_j,term_to_cellidscols_jt)

  matdata = (term_to_cellmat,term_to_cellidsrows, term_to_cellidscols)
end


function matdata_jacobian(op::TransientFEOperatorFromWeakForm,t::Real,uh::FEFunction,uh_t::FEFunction)
  Uh = evaluate(get_trial(op),nothing)
  du = get_cell_shapefuns_trial(Uh)
  v = get_cell_shapefuns(get_test(op))
  matdata = collect_cell_matrix(op.jac(t,uh,uh_t,du,v))
end

function matdata_jacobian_t(op::TransientFEOperatorFromWeakForm,t::Real,uh::FEFunction,uh_t::FEFunction,duht_du::Real)
  Uh = evaluate(get_trial(op),nothing)
  du_t = get_cell_shapefuns_trial(Uh)
  v = get_cell_shapefuns(get_test(op))
  matdata = collect_cell_matrix(duht_du*op.jac_t(t,uh,uh_t,du_t,v))
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
  jacobian_and_jacobian_t!(J,op,0.0,uh,uh,1.0,cache)
  @test isa(J,AbstractMatrix)
  cache = update_cache!(cache,op,0.0)
  true
end
