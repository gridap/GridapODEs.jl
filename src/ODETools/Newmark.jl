"""
Newmark-beta 2nd order ODE solver
"""
struct Newmark <: ODESolver
  nls::NonlinearSolver
  dt::Float64
  γ::Float64
  β::Float64
end

function solve_step!(
  u1::AbstractVector,
  solver::Newmark,
  op::ODEOperator,
  u0::AbstractVector,
  v0::AbstractVector,
  a0::AbstractVector,
  t0::Real,
  cache) # -> (uF,tF)
  
  dt = solver.dt
  t1 = t0+dt
  
  if cache === nothing
    ode_cache = allocate_cache(op)
    v1 = similar(v0)
    a1 = similar(a0)
    nl_cache = nothing
  else
    ode_cache, v1, a1, nl_cache = cache
  end
  
  ode_cache = update_cache!(ode_cache,op,t1)  
  nlop = NewmarkNonlinearOperator(op,t1,dt,u0,v0,a0,ode_cache,v1,a1)
  nl_cache = solve!(u1,solver.nls,nlop,nl_cache)  
  cache = (ode_cache, v1, a1, nl_cache)

  return (u1,t1,cache)
  
end

"""
Nonlinear operator that represents the Newmark nonlinear operator at a
given time step, i.e., A(t,u_n+1,v_n+1,a_n+1)
"""
struct ThetaMethodNonlinearOperator <: NonlinearOperator
  odeop::ODEOperator
  t1::Float64
  dt::Float64
  u0::AbstractVector
  v0::AbstractVector
  a0::AbstractVector
  ode_cache
  v1::AbstractVector
  a1::AbstractVector
end

function residual!(b::AbstractVector,op::NewmarkNonlinearOperator,x::AbstractVector)
  u1 = x
  v1 = op.v1
  a1 = op.a1
  a1 = 1.0/(op.β*op.dt^2)*(u1-op.u0) - 1.0/(op.β*op.dt)*op.v0 - (1-2*op.β)/(2*op.β)*op.a0
  v1 = op.γ/(op.β*op.dt)*(u1-op.u0) + (1-op.γ/op.β)*op.v0 + op.dt*(1-op.γ/(2*op.β))*op.a0
  residual!(b,op.odeop,op.t1,u1,v1,a1,op.ode_cache)
end

function jacobian!(A::AbstractMatrix,op::NewmarkNonlinearOperator,x::AbstractVector)
  u1 = x
  v1 = op.v1
  a1 = op.a1
  a1 = 1.0/(op.β*op.dt^2)*(u1-op.u0) - 1.0/(op.β*op.dt)*op.v0 - (1-2*op.β)/(2*op.β)*op.a0
  v1 = op.γ/(op.β*op.dt)*(u1-op.u0) + (1-op.γ/op.β)*op.v0 + op.dt*(1-op.γ/(2*op.β))*op.a0
  z = zero(eltype(A))
  fill_entries!(A,z)
  jacobian_and_jacobian_t!(A,op.odeop,op.t1,u1,v1,a1,op.γ/(op.β*op.dt)*(u1-op.u0),1.0/(op.β*op.dt^2),op.ode_cache)
end

function allocate_residual(op::NewmarkNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,x,op.ode_cache)
end

function allocate_jacobian(op::NewmarkNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,x,op.ode_cache)
end

function zero_initial_guess(op::NewmarkNonlinearOperator)
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end