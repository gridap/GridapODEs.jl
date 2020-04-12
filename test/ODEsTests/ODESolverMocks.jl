using Gridap.Algebra: residual
using Gridap.Algebra: jacobian
using Gridap.Algebra: fill_entries!
import Gridap.Algebra: NonlinearSolver
import Gridap.Algebra: NonlinearOperator
import Gridap.Algebra: solve!
import GridapTimeStepper.ODETools: solve_step!
import GridapTimeStepper.ODETools: ODESolver
import GridapTimeStepper.ODETools: zero_initial_guess
import GridapTimeStepper.ODETools: residual!
import GridapTimeStepper.ODETools: jacobian!
import GridapTimeStepper.ODETools: solve!
import GridapTimeStepper.ODETools: allocate_residual
import GridapTimeStepper.ODETools: allocate_jacobian

struct OperatorMock <: NonlinearOperator
  odeop
  tf::Float64
  dt::Float64
  u0::AbstractVector
  cache
end

function OperatorMock(odeop::ODEOperator,tf::Real,dt::Real,u0::AbstractVector)
  cache = nothing
  OperatorMock(odeop,tf,dt,u0,cache)
end

function residual!(b::AbstractVector,op::OperatorMock,x::AbstractVector)
  uf = x
  uf_t = (x-op.u0)/op.dt
  residual!(b,op.odeop,op.tf,uf,uf_t,op.cache)
end

function jacobian!(A::AbstractMatrix,op::OperatorMock,x::AbstractVector)
  uf = x
  uf_t = (x-op.u0)/op.dt
  fill_entries!(A,0.0)
  jacobian!(A,op.odeop,op.tf,uf,uf_t,op.cache)
  jacobian_t!(A,op.odeop,op.tf,uf,uf_t,(1/op.dt),op.cache)
end

function allocate_residual(op::OperatorMock,x::AbstractVector)
  allocate_residual(op.odeop,x,op.cache)
end

function allocate_jacobian(op::OperatorMock,x::AbstractVector)
  allocate_jacobian(op.odeop,x,op.cache)
end

function zero_initial_guess(::Type{T},op::OperatorMock) where T
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end

struct NLSolverMock <: NonlinearSolver
end

function solve!(x::AbstractVector,nls::NLSolverMock,nlop::NonlinearOperator)
  r = residual(nlop,x)
  J = jacobian(nlop,x)
  dx = inv(J)*(-r)
  x.= x.+dx
  cache = (r,J,dx)
end

function solve!(x::AbstractVector,nls::NLSolverMock,nlop::NonlinearOperator,cache)
  r, J, dx = cache
  residual!(r, nlop, x)
  jacobian!(J, nlop, x)
  dx = inv(J)*(-r)
  x.= x.+dx
end

struct ODESolverMock <: ODESolver
  nls::NLSolverMock
  dt::Float64
end

function solve_step!(
  uf::AbstractVector,solver::ODESolverMock,op::ODEOperator,u0::AbstractVector,t0::Real, ode_cache, nl_cache) # -> (uF,tF)


  dt = solver.dt
  tf = t0+dt
  update_cache!(ode_cache,op,tf)
  nlop = OperatorMock(op,tf,dt,u0,ode_cache)

  if (nl_cache==nothing)
    nl_cache = solve!(uf,solver.nls,nlop)
  else
    nl_cache = solve!(uf,solver.nls,nlop,nl_cache)
  end

  return (uf, tf, ode_cache, nl_cache)
end
