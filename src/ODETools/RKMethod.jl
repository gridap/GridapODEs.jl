struct ButcherTableType{Kind} end

"""
Butcher table
"""
struct ButcherTable
  s::Int # stages
  p::Int # embedded order
  q::Int # order
  a::Matrix # A_ij
  b::Vector # b_j
  c::Vector # c_i
  d::Vector # d_j (embedded)
  type::ButcherTableType # identifier
  function ButcherTable(type::ButcherTableType)
    createButcherTable(type)
  end
end

"""
Runge-Kutta ODE solver
"""
struct RKMethod <: ODESolver
  nls::NonlinearSolver
  dt::Float64
  bt::ButcherTable
  function RK(nls,dt,bt)
    new(nls,dt,bt)
  end
end

function createButcherTable(type::ButcherTableType{:SDIRK_2_1_2})
  s = 2
  p = 1
  q = 2
  a = [1.0 0.0; -1.0 1.0]
  b = [0.5 0.5]
  c = [1.0 0.0]
  d = [1.0 0.0]
  ButcherTable(s,p,q,a,b,c,d,type)
end

function solve_step!(uf::AbstractVector,
  solver::RKMethod,
  op::ODEOperator,
  u0::AbstractVector,
  t0::Real,
  cache)

  # Unpack variables
  dt = solver.dt
  s = solver.bt.s
  a = solver.bt.a
  b = solver.bt.b
  c = solver.bt.c
  d = solver.bt.d

  # Compute intermediate stages
  for i in 1:s

    # Create cache if not there
    if cache === nothing
      ode_cache = allocate_cache(op)
      vi = similar(u0)
      fi = [similar(u0)]
      nl_cache = nothing
    else
      ode_cache, vi, fi, nl_cache = cache
    end

    # allocate space to store the RHS at i
    if (length(fi) < i)
      push!(fi,similar(u0))
    end

    # solve at stage i
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    nlop = RKMethodNonlinearOperator(op,ti,dt,u0,ode_cache,vi,fi,i,a)
    nl_cache = solve!(uf,solver.nls,nlop,nl_cache)
    fi[i] = nlop.fi[i]

    cache = (ode_cache, vi, fi, nl_cache)

  end

  # update
  uf = u0
  for i in 1:s
    uf = uf + dt*b[i]*fi[i]
  end

  tf = t0 + dt
  return (uf,tf,cache)

end

"""
Nonlinear operator that represents the Runge-Kutta nonlinear operator at a
given time step, i.e., A(t,u_i,(u_i-u_n)/dt)
"""
struct RKMethodNonlinearOperator <: NonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vi::AbstractVector
  fi::AbstractVector
  i::Int
  a::Matrix
end

function residual!(b::AbstractVector,op::RKMethodNonlinearOperator,x::AbstractVector)
  # ∂u/∂t - a_ii * f(ui,ti) - ∑_{j<i} a_ij * f(uj,tj) = 0
  # Res_i = [1/a_ii * ∂u/∂t - f(ui,ti)]
  # Res_ij = - a_ij/a_ii * f(uj,ti)
  # Res_i + ∑_{j<i} Res_ij = 0
  ui = x
  vi = op.vi
  vi = (X-op.u0)/(op.a[op.i,op.i]*op.dt)
  residual!(b,op.odeop,op.ti,ui,vi,op.ode_cache)
  op.fi[op.i] = (vi-b)*op.a[op.i,op.i] # store fi for future stages
  for j in 1:op.i
    b = b - op.a[op.i,op.j]/op.a[op.i,op.i] * op.fi[op.j]
  end
end

function jacobian!(A::AbstractMatrix,op::RKMethodNonlinearOperator,x::AbstractVector)
  ui = x
  vi = op.vi
  vi = (X-op.u0)/(op.a[op.i,op.i]*op.dt)
  z = zero(eltype(A))
  fill_entries!(A,z)
  jacobian_and_jacobian_t!(A,op.odeop,op.tθ,uF,vi,(1/(op.a[op.i,op.i]*op.dt)),op.ode_cache)
end

function allocate_residual(op::RKMethodNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,x,op.ode_cache)
end

function allocate_jacobian(op::RKMethodNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,x,op.ode_cache)
end

function zero_initial_guess(op::RKMethodNonlinearOperator)
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end
