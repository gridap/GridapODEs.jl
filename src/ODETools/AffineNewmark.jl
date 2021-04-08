function solve_step!(
  x1::NTuple{3,AbstractVector},
  solver::Newmark,
  op::AffineODEOperator,
  x0::NTuple{3,AbstractVector},
  t0::Real,
  cache) # -> (uF,tF)

  dt = solver.dt
  γ = solver.γ
  β = solver.β
  t1 = t0+dt
  u0, v0, a0 = x0
  u1, v1, a1 = x1
  newmatrix = true

  if cache === nothing
    newmark_cache = allocate_cache(op,v0,a0)
    A, b = _allocate_matrix_and_vector(op,x0)
    affOp_cache = (A,b,)
  else
    newmark_cache, affOp_cache = cache
  end

  (v,a, ode_cache) = newmark_cache
  ode_cache = update_cache!(ode_cache,op,t1)
  newmark_affOp = NewmarkAffineOperator(op,t1,dt,γ,β,(u0,v0,a0),newmark_cache, affOp_cache)
  _matrix_and_vector!(newmark_affOp,u1)
  A,b,l_cache = newmark_affOp.affop_cache
  affOp = AffineOperator(A,b)
  l_cache = solve!(u1,solver.nls,affOp,l_cache,newmatrix)

  v1 = γ/(β*dt)*(u1-u0) + (1-γ/β)*v0 + dt*(1-γ/(2*β))*a0
  a1 = 1.0/(β*dt^2)*(u1-u0) - 1.0/(β*dt)*v0 - (1-2*β)/(2*β)*a0

  affOp_cache = A,b,l_cache
  cache = (newmark_cache, affOp_cache)
  x1 = (u1,v1,a1)

  return (x1,t1,cache)

end

"""
Affine operator that represents the Newmark Affine operator at a
given time step, i.e., M(t)(u_n+1-u_n)/dt + K(t)u_n+1 + b(t)
"""
struct NewmarkAffineOperator <: AffineOperator
  odeop::AffineODEOperator
  t1::Float64
  dt::Float64
  γ::Float64
  β::Float64
  x0::NTuple{3,AbstractVector}
  ode_cache
  affop_cache
end

function _matrix_and_vector!(affOp::NewmarkAffineOperator,x::AbstractVector)
  A, b, l_cache = affOp.affop_cache
  jacobian!(A,affOp,x)
  residual!(b,affOp,x)
  affOp.affop_cache = A, b, l_cache
  nothing
end

function residual!(b::AbstractVector,op::NewmarkAffineOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  residual!(b,op.odeop,op.t1,(u1,v1,a1),cache)
  b .*= -1.0
end

function jacobian!(A::AbstractMatrix,op::NewmarkAffineOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  z = zero(eltype(A))
  fill_entries!(A,z)
  jacobians!(A,op.odeop,op.t1,(u1,v1,a1),(1.0,op.γ/(op.β*op.dt),1.0/(op.β*op.dt^2)),cache)
end

function _mass_matrix!(A::AbstractMatrix,op::NewmarkAffineOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  z = zero(eltype(A))
  fill_entries!(A,z)
  jacobian!(A,op.odeop,op.t1,(u1,v1,a1),3,1.0/(op.β*op.dt^2),cache)
end

function allocate_residual(op::NewmarkAffineOperator,x::AbstractVector)
  v1, a1, cache = op.ode_cache
  allocate_residual(op.odeop,x,cache)
end

function allocate_jacobian(op::NewmarkAffineOperator,x::AbstractVector)
  v1, a1, cache = op.ode_cache
  allocate_jacobian(op.odeop,x,cache)
end

function _allocate_matrix_and_vector(op::NewmarkAffineOperator,x::AbstractVector)
  A = allocate_jacobian(op,x)
  b = allocate_residual(op,x)
  return A,b
end
