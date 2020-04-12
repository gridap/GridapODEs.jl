# const ∂t = time_derivative

# TransientFEOPerator cache
struct TransientFEOperatorCache
  Uh::FESpace
  Uht::FESpace
end

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

function allocate_residual(op::TransientFEOperator,uh)#,cache)
  @notimplemented
end

function allocate_jacobian(op::TransientFEOperator,uh,cache)
  @notimplemented
end

"""
Idem as `residual!` of `ODEOperator`
"""
function residual!(b::AbstractVector,op::TransientFEOperator,t,uh,uht,cache)
  @notimplemented
end

"""
Idem as `residual!` of `ODEOperator`
"""
function jacobian!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht,cache)
  @notimplemented
end

"""
Idem as `jacobian_t!` of `ODEOperator`
"""
function jacobian_t!(A::AbstractMatrix,op::TransientFEOperator,t,uh,uht,duht_du,cache)
  @notimplemented
end

"""
Returns a `ODEOperator` wrapper of the `TransientFEOperator` that can be
straightforwardly used with the `ODETools` module.
"""
function get_algebraic_operator(feop::TransientFEOperator)
  ODEOpFromFEOp(feop)
end

"""
Allocates the `cache` of a `TransientFEOperator`, which includes a `FESpace`
for the trial space for the unknown and its time derivative. It is used to
pre-allocate the arrays of fixed Dirichlet values once and overwrite them
at every time step (if transient trial spaces are being used due to variable
strong Dirichlet boundary conditions)
"""
function allocate_cache(feop::TransientFEOperator)
  cache = _allocate_cache(get_trial(feop))
  # assem = get_assembler(feop)
  cache
end

"""
Updates the cache, i.e., the strong Dirichlet values when dealing with
time-dependent Dirichlet data.
"""
function update_cache!(cache,feop::TransientFEOperator,t::Real)
  _update_cache!(cache,get_trial(feop),t)
end

"""
Returns the assembler, which is constant for all time steps for a given FE
operator.

Note: adaptive FE spaces involve to generate new FE spaces and
corresponding operators, due to the ummutable approach in `Gridap`
"""
get_assembler(feop::TransientFEOperator) = @notimplemented

# Internal functions

function _allocate_cache(fesp::FESpace)
  Uh = fesp
  Uht = HomogeneousTrialFESpace(fesp)
  TransientFEOperatorCache(Uh,Uht)
end

function _allocate_cache(fesp::TransientTrialFESpace)
  Uh = HomogeneousTrialFESpace(fesp.space)
  Uht = HomogeneousTrialFESpace(fesp.space)
  TransientFEOperatorCache(Uh,Uht)
end

_update_cache!(cache,::FESpace,t) = nothing

function _update_cache!(cache,tfesp::TransientTrialFESpace,t::Real)
  Uh = cache.Uh; Uht = cache.Uht
  TrialFESpace!(Uh,0.0)
  TrialFESpace!(Uh,tfesp.dirichlet_t(t))
  fun = tfesp.dirichlet_t
  fun_t = ∂t(fun)
  TrialFESpace!(Uht,fun_t(t))
  cache #Uh, Uht
end


# Specializations

"""
Transient FE operator that is defined by a set of `TransientFETerm` (or `FETerm`)
"""
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
  assem_t = SparseMatrixAssembler(test,trial(0.0))
  TransientFEOperatorFromTerms(trial,∂t(trial),test,assem_t,terms...)
end

get_test(op::TransientFEOperatorFromTerms) = op.test

get_trial(op::TransientFEOperatorFromTerms) = op.trial

get_trial(op::TransientFEOperatorFromTerms,t) = op.trial(t)

function allocate_residual(op::TransientFEOperatorFromTerms,uh)#,cache)#,assem)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  _, cellids =
  collect_cell_residual(0.0,uh,uh,v,op.terms)
  allocate_vector(op.assem_t,cellids)
end

function residual!(b::AbstractVector,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,cache)#,assem)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  v = get_cell_basis(op.test)
  cellvecs, cellids = collect_cell_residual(t,uh,uh_t,v,op.terms)
  assemble_vector!(b,op.assem_t,cellvecs,cellids)
  b
end

function allocate_jacobian(op::TransientFEOperatorFromTerms,uh,cache)
  Uh = cache.Uh;  Uht = cache.Uht
  @assert is_a_fe_function(uh)
  du = get_cell_basis(Uh)
  v = get_cell_basis(op.test)
  _, cellidsrows, cellidscols = collect_cell_jacobian(0.0,uh,uh,du,v,op.terms)
  allocate_matrix(op.assem_t, cellidsrows, cellidscols)
end

function jacobian!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,cache)#,assem)
  Uh = cache.Uh;  Uht = cache.Uht
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  du = get_cell_basis(Uh)
  v = get_cell_basis(op.test)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian(t,uh,uh_t,du,v,op.terms)
  assemble_matrix_add!(A,op.assem_t, cellmats, cellidsrows, cellidscols)
  A
end

function jacobian_t!(A::AbstractMatrix,op::TransientFEOperatorFromTerms,
  t::Real,uh,uh_t,duht_du::Real,cache)#,assem)
  Uh = cache.Uh;  Uht = cache.Uht
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  du_t = get_cell_basis(Uh)
  v = get_cell_basis(op.test)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian_t(t,uh,uh_t,du_t,v,duht_du,op.terms)
  assemble_matrix_add!(A,op.assem_t, cellmats, cellidsrows, cellidscols)
  A
end

# Tester

function test_transient_fe_operator(op::TransientFEOperator,uh)
  odeop = get_algebraic_operator(op)
  @test isa(odeop,ODEOperator)
  cache = allocate_cache(odeop)
  V = get_test(op)
  @test isa(V,FESpace)
  U = get_trial(op)
  @test isa(U,Union{FESpace,TransientTrialFESpace})
  U0 = U(0.0)
  @test isa(U0,FESpace)
  r = allocate_residual(op,uh)
  @test isa(r,AbstractVector)
  residual!(r,op,0.0,uh,uh,cache)
  @test isa(r,AbstractVector)
  J = allocate_jacobian(op,uh,cache)
  @test isa(J,AbstractMatrix)
  jacobian!(J,op,0.0,uh,uh,cache)
  @test isa(J,AbstractMatrix)
  jacobian_t!(J,op,0.0,uh,uh,1.0,cache)
  @test isa(J,AbstractMatrix)
  true
end
