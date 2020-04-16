"""
θ-method ODE solver
"""
struct ThetaMethod <: ODESolver
  nls::NonlinearSolver
  dt::Float64
  θ::Float64
end

function solve_step!(uf::AbstractVector,
                     solver::ThetaMethod,
                     op::ODEOperator,
                     u0::AbstractVector,
                     t0::Real,
                     ode_cache,
                     nl_cache) # -> (uF,tF)

  # Build the nonlinear problem to solve at this step
  dt = solver.dt
  tf = t0+dt
  solver.θ == 0.0 ? dtθ = dt : dtθ = dt*solver.θ
  tθ = t0+dtθ
  # if 0.0 < solver.θ < 1.0
  #   v0 = (1-solver.θ)*u0
  # else
  #   v0 = u0
  # end
  update_cache!(ode_cache,op,tf)

  nlop = ThetaMethodNonlinearOperator(op,tθ,dtθ,u0,ode_cache)

  # Solve the nonlinear problem
  if (nl_cache==nothing)
    nl_cache = solve!(uf,solver.nls,nlop)
  else
    # solve!(uf,solver.nls,nlop,nl_cache)
    # @santiagobadia: What I am doing wrong here?
    # could we create a function with methods dispatching based on solver
    # linear or nonlinear? 
    x = copy(nlop.u0)
    b = allocate_residual(nlop,x)
    residual!(b,nlop,x)
    J = allocate_jacobian(nlop,x)
    jacobian!(J,nlop,x)
    uf = u0-J\b
  end

  if 0.0 < solver.θ < 1.0
    uf = uf*(1.0/solver.θ)-u0*((1-solver.θ)/solver.θ)
  end

  return (uf,tf,ode_cache,nl_cache)

end

# Struct representing the nonlinear algebraic problem to be solved at a given step
"""
Nonlinear operator that represents the θ-method nonlinear operator at a
given time step, i.e., A(t,u_n+θ,(u_n+θ-u_n)/dt)
"""
struct ThetaMethodNonlinearOperator <: NonlinearOperator
  odeop::ODEOperator
  tθ::Float64
  dtθ::Float64
  u0::AbstractVector
  ode_cache
end

function residual!(b::AbstractVector,op::ThetaMethodNonlinearOperator,x::AbstractVector)
  uθ = x
  vθ = (x-op.u0)/op.dtθ
  residual!(b,op.odeop,op.tθ,uθ,vθ,op.ode_cache)
end

function jacobian!(A::AbstractMatrix,op::ThetaMethodNonlinearOperator,x::AbstractVector)
  uF = x
  vF = (x-op.u0)/op.dtθ
  z = zero(eltype(A))
  fill_entries!(A,z)
  jacobian!(A,op.odeop,op.tθ,uF,vF,op.ode_cache)
  jacobian_t!(A,op.odeop,op.tθ,uF,vF,(1/op.dtθ),op.ode_cache)
end

function allocate_residual(op::ThetaMethodNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,x,op.ode_cache)
end

function allocate_jacobian(op::ThetaMethodNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,x,op.ode_cache)
end

function zero_initial_guess(::Type{T},op::ThetaMethodNonlinearOperator) where T
  x0 = similar(op.u0,T)
  fill!(x0,zero(T))
  x0
end
