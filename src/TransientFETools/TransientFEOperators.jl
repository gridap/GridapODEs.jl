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
function residual!(
  b::AbstractVector,
  op::TransientFEOperator,
  t::Real,
  xh::Union{AbstractVector,Tuple{Vararg{AbstractVector}}},
  cache)
  @abstractmethod
end

"""
Idem as `jacobian!` of `ODEOperator`
"""
function jacobian!(
  A::AbstractMatrix,
  op::TransientFEOperator,
  t::Real,
  xh::Union{AbstractVector,Tuple{Vararg{AbstractVector}}},
  i::Int,
  γᵢ::Real,
  cache)
  @abstractmethod
end

"""
Idem as `jacobians!` of `ODEOperator`
"""
function jacobians!(
  A::AbstractMatrix,
  op::TransientFEOperator,
  t::Real,
  x::Tuple{Vararg{AbstractVector}},
  γ::Tuple{Vararg{Real}},
  cache)
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
  jacs::Tuple{Vararg{Function}}
  assem_t::Assembler
  trials::Tuple{Vararg{Any}}
  test::FESpace
  order::Integer
end

function TransientConstantFEOperator(m::Function,a::Function,b::Function,
  trial,test)
  res(t,(u,ut),v) = -1.0 * b(v)
  jac(t,(u,ut),du,v) = a(du,v)
  jac_t(t,(u,ut),dut,v) = m(dut,v)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{Constant}(res,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end

function TransientConstantMatrixFEOperator(m::Function,a::Function,b::Function,
  trial,test)
  res(t,(u,ut),v) = m(ut,v) + a(u,v) - b(t,v)
  jac(t,(u,ut),du,v) = a(du,v)
  jac_t(t,(u,ut),dut,v) = m(dut,v)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{ConstantMatrix}(res,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end

function TransientAffineFEOperator(m::Function,a::Function,b::Function,
  trial,test)
  res(t,(u,ut),v) = m(t,ut,v) + a(t,u,v) - b(t,v)
  jac(t,(u,ut),du,v) = a(t,du,v)
  jac_t(t,(u,ut),dut,v) = m(t,dut,v)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{Affine}(res,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end

function TransientFEOperator(res::Function,jac::Function,jac_t::Function,
  trial,test)
  assem_t = SparseMatrixAssembler(trial,test)
  TransientFEOperatorFromWeakForm{Nonlinear}(res,(jac,jac_t),assem_t,(trial,∂t(trial)),test,1)
end


function TransientConstantFEOperator(m::Function,c::Function,a::Function,b::Function,
  trial,test)
  res(t,(u,ut,utt),v) = -1.0 * b(v)
  jac(t,(u,ut,utt),du,v) = a(du,v)
  jac_t(t,(u,ut,utt),dut,v) = c(dut,v)
  jac_tt(t,(u,ut,utt),dutt,v) = m(dutt,v)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  TransientFEOperatorFromWeakForm{Constant}(
    res,(jac,jac_t,jac_tt),assem_t,(trial,trial_t,trial_tt),test,2)
end

function TransientConstantMatrixFEOperator(m::Function,c::Function,a::Function,b::Function,
  trial,test)
  res(t,(u,ut,utt),v) = m(utt,v) + c(ut,v) + a(u,v) - b(t,v)
  jac(t,(u,ut,utt),du,v) = a(du,v)
  jac_t(t,(u,ut,utt),dut,v) = c(dut,v)
  jac_tt(t,(u,ut,utt),dutt,v) = m(dutt,v)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  TransientFEOperatorFromWeakForm{ConstantMatrix}(
    res,(jac,jac_t,jac_tt),assem_t,(trial,trial_t,trial_tt),test,2)
end

function TransientAffineFEOperator(m::Function,c::Function,a::Function,b::Function,
  trial,test)
  res(t,(u,ut,utt),v) = m(t,utt,v) + c(t,ut,v) + a(t,u,v) - b(t,v)
  jac(t,(u,ut,utt),du,v) = a(t,du,v)
  jac_t(t,(u,ut,utt),dut,v) = c(t,dut,v)
  jac_tt(t,(u,ut,utt),dutt,v) = m(t,dutt,v)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  TransientFEOperatorFromWeakForm{Affine}(
    res,(jac,jac_t,jac_tt),assem_t,(trial,trial_t,trial_tt),test,2)
end

function TransientFEOperator(res::Function,jac::Function,jac_t::Function,
  jac_tt::Function,trial,test)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  TransientFEOperatorFromWeakForm{Nonlinear}(
    res,(jac,jac_t,jac_tt),assem_t,(trial,trial_t,trial_tt),test,2)
end

function TransientFEOperator(res::Function,trial,test;order::Integer=1)
  jacs = ()
  for i in 1:order+1
    function jac_i(t,x,dxi,dv)
      function res_i(y)
        x[i] = y
        res(t,x,dv)
      end
      jacobian(res_i,x[i])
    end
    jacs = (jacs...,jac_i)
  end
  # jac(t,u,ut,du,dv) = jacobian(x->res(t,x,ut,dv),u)
  # jac_t(t,u,ut,dut,dv) = jacobian(xt->res(t,u,xt,dv),ut)
  TransientFEOperator(res,jacs,trial,test)
end

function SparseMatrixAssembler(
  trial::Union{TransientTrialFESpace,TransientMultiFieldTrialFESpace},
  test::FESpace)
  SparseMatrixAssembler(evaluate(trial,nothing),test)
end

get_assembler(op::TransientFEOperatorFromWeakForm) = op.assem_t

get_test(op::TransientFEOperatorFromWeakForm) = op.test

get_trial(op::TransientFEOperatorFromWeakForm) = op.trials[1]

get_order(op::TransientFEOperatorFromWeakForm) = op.order

function allocate_residual(op::TransientFEOperatorFromWeakForm,uh::FEFunction,cache)
  V = get_test(op)
  v = get_cell_shapefuns(V)
  xh = ()
  for i in 1:get_order(op)+1
    xh = (xh...,uh)
  end
  vecdata = collect_cell_vector(V,op.res(0.0,xh,v))
  allocate_vector(op.assem_t,vecdata)
end

function residual!(
  b::AbstractVector,
  op::TransientFEOperatorFromWeakForm,
  t::Real,
  xh::Tuple{Vararg{FEFunction}},
  cache)
  V = get_test(op)
  v = get_cell_shapefuns(V)
  vecdata = collect_cell_vector(V,op.res(t,xh,v))
  assemble_vector!(b,op.assem_t,vecdata)
  b
end

function allocate_jacobian(op::TransientFEOperatorFromWeakForm,uh::FEFunction,cache)
  xh = ()
  for i in 1:get_order(op)+1
    xh = (xh...,uh)
  end
  _matdata = ()
  for i in 1:get_order(op)+1
    _matdata = (_matdata...,matdata_jacobian(op,0.0,xh,i,0.0))
  end
  matdata = vcat_matdata(_matdata)
  allocate_matrix(op.assem_t,matdata)
end

function jacobian!(
  A::AbstractMatrix,
  op::TransientFEOperatorFromWeakForm,
  t::Real,
  xh::Tuple{Vararg{FEFunction}},
  i::Integer,
  γᵢ::Real,
  cache)
  matdata = matdata_jacobian(op,t,xh,i,γᵢ)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function jacobians!(
  A::AbstractMatrix,
  op::TransientFEOperatorFromWeakForm,
  t::Real,
  xh::Tuple{Vararg{FEFunction}},
  γ::Tuple{Vararg{Real}},
  cache)
  _matdata = ()
  for i in 1:get_order(op)+1
    _matdata = (_matdata...,matdata_jacobian(op,t,xh,i,γ[i]))
  end
  matdata = vcat_matdata(_matdata)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function vcat_matdata(_matdata)
  term_to_cellmat_j = ()
  term_to_cellidsrows_j = ()
  term_to_cellidscols_j = ()
  for j in 1:length(_matdata)
    term_to_cellmat_j = (term_to_cellmat_j...,_matdata[j][1])
    term_to_cellidsrows_j = (term_to_cellidsrows_j...,_matdata[j][2])
    term_to_cellidscols_j = (term_to_cellidscols_j...,_matdata[j][3])
  end

  term_to_cellmat = vcat(term_to_cellmat_j...)
  term_to_cellidsrows = vcat(term_to_cellidsrows_j...)
  term_to_cellidscols = vcat(term_to_cellidscols_j...)

  matdata = (term_to_cellmat,term_to_cellidsrows, term_to_cellidscols)
end


function matdata_jacobian(
  op::TransientFEOperatorFromWeakForm,
  t::Real,
  xh::Tuple{Vararg{FEFunction}},
  i::Integer,
  γᵢ::Real)
  Uh = evaluate(get_trial(op),nothing)
  V = get_test(op)
  du = get_cell_shapefuns_trial(Uh)
  v = get_cell_shapefuns(V)
  matdata = collect_cell_matrix(Uh,V,γᵢ*op.jacs[i](t,xh,du,v))
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
  residual!(r,op,0.0,(uh,uh),cache)
  @test isa(r,AbstractVector)
  J = allocate_jacobian(op,uh,cache)
  @test isa(J,AbstractMatrix)
  jacobian!(J,op,0.0,(uh,uh),1,1.0,cache)
  @test isa(J,AbstractMatrix)
  jacobian!(J,op,0.0,(uh,uh),2,1.0,cache)
  @test isa(J,AbstractMatrix)
  jacobians!(J,op,0.0,(uh,uh),(1.0,1.0),cache)
  @test isa(J,AbstractMatrix)
  cache = update_cache!(cache,op,0.0)
  true
end
