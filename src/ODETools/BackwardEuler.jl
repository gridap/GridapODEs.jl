"""
Backward Euler ODE solver
"""
struct BackwardEuler <: ODESolver
  nls::NonlinearSolver
  dt::Float64
end

function solve_step!(
  uf::AbstractVector,solver::BackwardEuler,op::ODEOperator,u0::AbstractVector,t0::Real,cache) # -> (uF,tF,cache)

  if cache === nothing
    ode_cache = allocate_cache(op)
    vF = similar(u0)
    nl_cache = nothing
  else
    ode_cache, vF, nl_cache = cache
  end

  # Build the nonlinear problem to solve at this step
  dt = solver.dt
  tf = t0+dt
  ode_cache = update_cache!(ode_cache,op,tf)
  nlop = BackwardEulerNonlinearOperator(op,tf,dt,u0,ode_cache,vF)

  # Solve the nonlinear problem
  if (nl_cache==nothing)
    nl_cache = solve!(uf,solver.nls,nlop)
  else
    solve!(uf,solver.nls,nlop,nl_cache)
  end

  cache = (ode_cache, vF, nl_cache)

  return (uf,tf,cache)
end

# Struct representing the nonlinear algebraic problem to be solved at a given step
"""
Nonlinear operator that represents the Backward Euler nonlinear operator at a
given time step, i.e., A(t,u_n+1,(u_n+1-u_n)/dt)
"""
struct BackwardEulerNonlinearOperator <: NonlinearOperator
  odeop::ODEOperator
  tF::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vF:AbstractVector
end

function residual!(b::AbstractVector,op::BackwardEulerNonlinearOperator,x::AbstractVector)
  uF = x
  vF = op.vF
  vF .= (x .- op.u0) ./ op.dt
  residual!(b,op.odeop,op.tF,uF,vF,op.ode_cache)
end

function jacobian!(A::AbstractMatrix,op::BackwardEulerNonlinearOperator,x::AbstractVector)
  uF = x
  vF = op.vF
  vF .= (x .- op.u0) ./ op.dt
  z = zero(eltype(A))
  fill_entries!(A,z)
  jacobian!(A,op.odeop,op.tF,uF,vF,op.ode_cache)
  jacobian_t!(A,op.odeop,op.tF,uF,vF,(1/op.dt),op.ode_cache)
end

function allocate_residual(op::BackwardEulerNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,x,op.ode_cache)
end

function allocate_jacobian(op::BackwardEulerNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,x,op.ode_cache)
end

function zero_initial_guess(::Type{T},op::BackwardEulerNonlinearOperator) where T
  x0 = similar(op.u0,T)
  fill!(x0,zero(T))
  x0
end
