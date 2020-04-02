import Gridap.Algebra: NonLinearSolver
import Gridap.Algebra: NonLinearOperator
import Gridap.Algebra: solve!
import GridapTimeStepper.ODETools: solve_step!
import GridapTimeStepper.ODETools: allocate_cache
import GridapTimeStepper.ODETools: ODESolver
using Gridap.Algebra: residual
using Gridap.Algebra: jacobian
import GridapTimeStepper.ODETools: zero_initial_guess
import GridapTimeStepper.ODETools: residual!
import GridapTimeStepper.ODETools: jacobian!
import GridapTimeStepper.ODETools: solve!
# using GridapTimeStepper.ODETools: solve!
# import Gridap.Algebra: solve!

function fill_entries!(J,v)
  J .= zero(eltype(J))
end

struct OperatorMock <: NonLinearOperator
  odeop
  tf::Float64
  dt::Float64
  u0::AbstractVector
end

function residual!(b::AbstractVector,op::OperatorMock,x::AbstractVector)
  uf = x
  uf_t = (x-op.u0)/op.dt
  residual!(b,op.odeop,op.tf,uf,uf_t)
end

function jacobian!(A::AbstractMatrix,op::OperatorMock,x::AbstractVector)
  uf = x
  uf_t = (x-op.u0)/op.dt
  fill_entries!(A,zero(eltype(A)))
  jacobian_unknown!(A,op.odeop,op.tf,uf,uf_t)
  jacobian_unknown_t!(A,op.odeop,op.tf,uf,uf_t,(1/op.dt))
end

function allocate_residual(op::OperatorMock,x::AbstractVector)
  allocate_residual(op.odeop,x,x)
end

function allocate_jacobian(op::OperatorMock,x::AbstractVector)
  allocate_jacobian(op.odeop,x,x)
end

function zero_initial_guess(::Type{T},op::OperatorMock) where T
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end

struct NLSolverMock <: NonLinearSolver
end

function solve!(x::AbstractVector,nls::NLSolverMock,nlop::NonLinearOperator)
  r = residual(nlop,x)
  J = jacobian(nlop,x)
  dx = inv(J)*(-r)
  x.= x.+dx
  cache = (r,J,dx)
end

function solve!(x::AbstractVector,nls::NLSolverMock,nlop::NonLinearOperator,cache)
  r, J, dx = cache
  residual!(r, op, x)
  jacobian!(J, op, x)
  dx = inv(J)*(-r)
  x += dx
end

struct ODESolverMock <: ODESolver
  nls::NLSolverMock
  dt::Float64
end

function solve_step!(
  uf::AbstractVector,solver::ODESolverMock,op::ODEOperator,u0::AbstractVector,t0::Real, cache) # -> (uF,tF)

  dt = solver.dt
  tf = t0+dt
  nlop = OperatorMock(op,tf,dt,u0)

  if (cache==nothing)
    cache = solve!(uf,solver.nls,nlop)
  else
    cache = solve!(uf,solver.nls,nlop,cache)
  end

  return (uf, tf, cache)
end

function allocate_cache(
  solver::ODESolverMock,op::ODEOperator,u0::AbstractVector,t0::Real)
  r = allocate_residual(solver.nls,x)
  J = allocate_jacobian(solver.nls,x)
  dx = copy(u0)
  (r, J, dx)

end
