import Gridap.Algebra: NonLinearSolver
import Gridap.Algebra: NonLinearOperator
import Gridap.Algebra: solve!
import GridapTimeStepper.ODETools: solve_step!
import GridapTimeStepper.ODETools: allocate_cache
import GridapTimeStepper.ODETools: ODESolver

struct NLSolverMock <: NonLinearSolver end

  # uF, cache = solve!(uF,solver.nls,nlop,cache) # TODO reuse the cache
solve!(uF::AbstractVector,nls::NLSolverMock,nlop::NonLinearOperator,cache) = (nlop.u0+nlop.dt*nlop.u0,cache)

struct OperatorMock <: NonLinearOperator
  op
  tF::Float64
  dt::Float64
  u0::AbstractVector
end

struct ODESolverMock <: ODESolver
  nls::NLSolverMock
  dt::Float64
end

function solve_step!(
  uF::AbstractVector,solver::ODESolverMock,op::ODEOperator,u0::AbstractVector,t0::Real, cache) # -> (uF,tF)

  dt = solver.dt
  tF = t0+dt
  nlop = OperatorMock(op,tF,dt,u0)

  uF, cache = solve!(uF,solver.nls,nlop,cache) # TODO reuse the cache

  return (uF, tF, cache)
end

function allocate_cache(
  solver::ODESolverMock,op::ODEOperator,u0::AbstractVector,t0::Real)
  nothing
end
