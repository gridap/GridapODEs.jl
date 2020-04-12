# Now, we need an abstract type representing a numerical discretization scheme
# for the ODE
abstract type ODESolver <: GridapType end

function solve_step!(
  uF::AbstractVector,solver::ODESolver,op::ODEOperator,u0::AbstractVector,t0::Real,op_cache,cache) # -> (uF,tF,cache)
  @abstractmethod
end

function solve(
  solver::ODESolver,op::ODEOperator,u0::AbstractVector,t0::Real,tf::Real)
  GenericODESolution(solver,op,u0,t0,tf)
end

function test_ode_solver(solver::ODESolver,op::ODEOperator,u0,t0,tf)
  solution = solve(solver,op,u0,t0,tf)
  test_ode_solution(solution)
end

# I think we can create the machinery for ODE solvers easily, using more or less
# what we have above. But there is still a layer that combines it with the FE
# machinery in Gridap.

include("BackwardEuler.jl")
