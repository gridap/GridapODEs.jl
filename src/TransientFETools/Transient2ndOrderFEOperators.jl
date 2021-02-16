"""
Transient FE operator that is defined by a 2nd order transient Weak form
"""
struct Transient2ndOrderFEOperatorFromWeakForm{C} <: TransientFEOperator{C}
  res::Function
  jac::Function
  jac_t::Function
  jac_tt::Function
  assem_t::Assembler
  trial
  trial_t
  trial_tt
  test::FESpace
end

function TransientConstantFEOperator(m::Function,c::Function,a::Function,b::Function,
  trial,test)
  res(t,u,ut,utt,v) = m(utt,v) + c(ut,v) + a(u,v) - b(v)
  jac(t,u,ut,utt,du,v) = a(du,v)
  jac_t(t,u,ut,utt,dut,v) = c(dut,v)
  jac_tt(t,u,ut,utt,dutt,v) = m(dutt,v)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  Transient2ndOrderFEOperatorFromWeakForm{Constant}(
    res,jac,jac_t,jac_tt,assem_t,trial,trial_t,trial_tt,test)
end

function TransientAffineFEOperator(m::Function,c::Function,a::Function,b::Function,
  trial,test)
  res(t,u,ut,utt,v) = m(t,utt,v) + c(t,ut,v) + a(t,u,v) - b(t,v)
  jac(t,u,ut,utt,du,v) = a(t,du,v)
  jac_t(t,u,ut,utt,dut,v) = c(t,dut,v)
  jac_tt(t,u,ut,utt,dutt,v) = m(t,dutt,v)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  Transient2ndOrderFEOperatorFromWeakForm{Affine}(
    res,jac,jac_t,jac_tt,assem_t,trial,trial_t,trial_tt,test)
end

function TransientFEOperator(res::Function,jac::Function,jac_t::Function,
  jac_tt::Function,trial,test)
  assem_t = SparseMatrixAssembler(trial,test)
  trial_t = ∂t(trial)
  trial_tt = ∂t(trial_t)
  Transient2ndOrderFEOperatorFromWeakForm{Nonlinear}(
    res,jac,jac_t,jac_tt,assem_t,trial,trial_t,trial_tt,test)
end
function Transient2ndOrderFEOperator(res::Function,trial,test)
  jac(t,u,ut,utt,du,dv) = jacobian(x->res(t,x,ut,utt,dv),u)
  jac_t(t,u,ut,utt,dut,dv) = jacobian(xt->res(t,u,xt,utt,dv),ut)
  jac_tt(t,u,ut,utt,dutt,dv) = jacobian(xtt->res(t,u,ut,xtt,dv),utt)
  TransientFEOperator(res,jac,jac_t,jac_tt,trial,test)
end

get_assembler(op::Transient2ndOrderFEOperatorFromWeakForm) = op.assem_t

get_test(op::Transient2ndOrderFEOperatorFromWeakForm) = op.test

get_trial(op::Transient2ndOrderFEOperatorFromWeakForm) = op.trial

function allocate_residual(op::Transient2ndOrderFEOperatorFromWeakForm,uh::FEFunction,cache)
  v = get_cell_shapefuns(get_test(op))
  vecdata = collect_cell_vector(op.res(0.0,uh,uh,uh,v))
  allocate_vector(op.assem_t,vecdata)
end

function residual!(b::AbstractVector,op::Transient2ndOrderFEOperatorFromWeakForm,
  t::Real,uh::FEFunction,uh_t::FEFunction,uh_tt::FEFunction,cache)
  v = get_cell_shapefuns(get_test(op))
  vecdata = collect_cell_vector(op.res(t,uh,uh_t,uh_tt,v))
  assemble_vector!(b,op.assem_t,vecdata)
  b
end

function allocate_jacobian(op::Transient2ndOrderFEOperatorFromWeakForm,uh::FEFunction,cache)
  matdata_j = matdata_jacobian(op,0.0,uh,uh,uh)
  matdata_jt = matdata_jacobian_t(op,0.0,uh,uh,uh,0.0)
  matdata_jtt = matdata_jacobian_tt(op,0.0,uh,uh,uh,0.0)
  matdata = vcat_matdata(matdata_j,matdata_jt,matdata_jtt)
  allocate_matrix(op.assem_t,matdata)
end

function jacobian!(A::AbstractMatrix,op::Transient2ndOrderFEOperatorFromWeakForm,
  t::Real,uh::FEFunction,uh_t::FEFunction,uh_tt::FEFunction,cache)
  matdata = matdata_jacobian(op,t,uh,uh_t,uh_tt)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function jacobian_t!(A::AbstractMatrix,op::Transient2ndOrderFEOperatorFromWeakForm,
  t::Real,uh,uh_t,uh_tt,duht_du::Real,cache)
  matdata = matdata_jacobian_t(op,t,uh,uh_t,uh_tt,duht_du)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function jacobian_tt!(A::AbstractMatrix,op::Transient2ndOrderFEOperatorFromWeakForm,
  t::Real,uh,uh_t,uh_tt,duhtt_du::Real,cache)
  matdata = matdata_jacobian_tt(op,t,uh,uh_t,uh_tt,duhtt_du)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function jacobian_and_jacobian_t!(A::AbstractMatrix,op::Transient2ndOrderFEOperatorFromWeakForm,
  t::Real,uh,uh_t,uh_tt,duht_du::Real,duhtt_du::Real,cache)
  matdata_j = matdata_jacobian(op,t,uh,uh_t,uh_tt)
  matdata_jt = matdata_jacobian_t(op,t,uh,uh_t,uh_tt,duht_du)
  matdata_jtt = matdata_jacobian_tt(op,t,uh,uh_t,uh_tt,duhtt_du)
  matdata = vcat_matdata(matdata_j,matdata_jt,matdata_jtt)
  assemble_matrix_add!(A,op.assem_t, matdata)
  A
end

function vcat_matdata(matdata_j,matdata_jt,matdata_jtt)
  term_to_cellmat_j, term_to_cellidsrows_j, term_to_cellidscols_j = matdata_j
  term_to_cellmat_jt, term_to_cellidsrows_jt, term_to_cellidscols_jt = matdata_jt
  term_to_cellmat_jtt, term_to_cellidsrows_jtt, term_to_cellidscols_jtt = matdata_jtt

  term_to_cellmat = vcat(term_to_cellmat_j,term_to_cellmat_jt,term_to_cellmat_jtt)
  term_to_cellidsrows = vcat(term_to_cellidsrows_j,term_to_cellidsrows_jt,term_to_cellidsrows_jtt)
  term_to_cellidscols = vcat(term_to_cellidscols_j,term_to_cellidscols_jt,term_to_cellidscols_jtt)

  matdata = (term_to_cellmat,term_to_cellidsrows, term_to_cellidscols)
end


function matdata_jacobian(op::Transient2ndOrderFEOperatorFromWeakForm,t::Real,
  uh::FEFunction,uh_t::FEFunction,uh_tt::FEFunction)
  Uh = evaluate(get_trial(op),nothing)
  du = get_cell_shapefuns_trial(Uh)
  v = get_cell_shapefuns(get_test(op))
  matdata = collect_cell_matrix(op.jac(t,uh,uh_t,uh_tt,du,v))
end

function matdata_jacobian_t(op::Transient2ndOrderFEOperatorFromWeakForm,t::Real,
  uh::FEFunction,uh_t::FEFunction,uh_tt::FEFunction,duht_du::Real)
  Uh = evaluate(get_trial(op),nothing)
  du_t = get_cell_shapefuns_trial(Uh)
  v = get_cell_shapefuns(get_test(op))
  matdata = collect_cell_matrix(duht_du*op.jac_t(t,uh,uh_t,uh_tt,du_t,v))
end

function matdata_jacobian_tt(op::Transient2ndOrderFEOperatorFromWeakForm,t::Real,
  uh::FEFunction,uh_t::FEFunction,uh_tt::FEFunction,duhtt_du::Real)
  Uh = evaluate(get_trial(op),nothing)
  du_tt = get_cell_shapefuns_trial(Uh)
  v = get_cell_shapefuns(get_test(op))
  matdata = collect_cell_matrix(duhtt_du*op.jac_tt(t,uh,uh_t,uh_tt,du_tt,v))
end

"""
Returns a `SecondOrderODEOperator` wrapper of the `Transient2ndOrderFEOperator` that can be
straightforwardly used with the `ODETools` module.
"""
function get_algebraic_operator(feop::Transient2ndOrderFEOperatorFromWeakForm{C}) where C
  SecondOrderODEOpFromFEOp{C}(feop)
end

# Tester

function test_transient_2ndOrder_fe_operator(op::TransientFEOperator,uh)
  odeop = get_algebraic_operator(op)
  @test isa(odeop,SecondOrderODEOperator)
  cache = allocate_cache(op)
  V = get_test(op)
  @test isa(V,FESpace)
  U = get_trial(op)
  U0 = U(0.0)
  @test isa(U0,FESpace)
  r = allocate_residual(op,uh,cache)
  @test isa(r,AbstractVector)
  residual!(r,op,0.0,uh,uh,uh,cache)
  @test isa(r,AbstractVector)
  J = allocate_jacobian(op,uh,cache)
  @test isa(J,AbstractMatrix)
  jacobian!(J,op,0.0,uh,uh,uh,cache)
  @test isa(J,AbstractMatrix)
  jacobian_t!(J,op,0.0,uh,uh,uh,1.0,cache)
  @test isa(J,AbstractMatrix)
  jacobian_tt!(J,op,0.0,uh,uh,uh,1.0,cache)
  @test isa(J,AbstractMatrix)
  jacobian_and_jacobian_t!(J,op,0.0,uh,uh,uh,1.0,1.0,cache)
  @test isa(J,AbstractMatrix)
  cache = update_cache!(cache,op,0.0)
  true
end
