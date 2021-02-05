# Now, we need an abstract type representing a numerical discretization scheme
# for the ODE
"""
Represents a map that given (t_n,u_n) returns (t_n+1,u_n+1) and cache for the
corresponding `ODEOperator` and `NonlinearOperator`
"""
abstract type ODESolver <: GridapType end

function solve_step!(
  uF::AbstractVector,solver::ODESolver,op::ODEOperator,u0::AbstractVector,t0::Real,cache) # -> (uF,tF,cache)
  @abstractmethod
end

# Default API

function solve_step!(
  uF::AbstractVector,solver::ODESolver,op::ODEOperator,u0::AbstractVector,t0::Real) # -> (uF,tF,cache)
  solve_step!(uF,solver,op,u0,t0,nothing)
end

function solve(
  solver::ODESolver,op::ODEOperator,u0::AbstractVector,t0::Real,tf::Real)
  GenericODESolution(solver,op,u0,t0,tf)
end

# Default API 2nd Order ODE

function solve_step!(
  uF::AbstractVector,
  solver::ODESolver,
  op::ODEOperator,
  u0::AbstractVector,
  v0::AbstractVector,
  a0::AbstractVector,
  t0::Real) # -> (uF,tF,cache)
  solve_step!(uF,solver,op,u0,v0,a0,t0,nothing)
end

function solve(
  solver::ODESolver,
  op::ODEOperator,
  u0::AbstractVector,
  v0::AbstractVector,
  a0::AbstractVector,
  t0::Real,
  tf::Real)
  Generic2ndOrderODESolution(solver,op,u0,v0,a0,t0,tf)
end

# testers

function test_ode_solver(solver::ODESolver,op::ODEOperator,u0,t0,tf)
  solution = solve(solver,op,u0,t0,tf)
  test_ode_solution(solution)
end

# Specialization

include("ForwardEuler.jl")

include("ThetaMethod.jl")

include("AffineThetaMethod.jl")

include("RungeKutta.jl")

include("Newmark.jl")
