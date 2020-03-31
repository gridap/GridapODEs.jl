# Now, we need an abstract type representing a numerical discretization scheme for the ODE
#
abstract type ODESolver <: GridapType end

# The solver is defined by how a single step is solved
function solve_step!(
  uF::AbstractVector,solver::ODESolver,op::ODEOperator,u0::AbstractVector,t0::Real,cache) # -> (uF,tF)
  @abstractmethod
end

function solve_step!(
  uF::AbstractVector,solver::ODESolver,op::ODEOperator,u0::AbstractVector,t0::Real) # -> (uF,tF),cache
  @abstractmethod
end

# @santiagobadia : Agreed

# Assuming that solve_step! is defined, we can define the following solve function for all ODESolver objects:
# Given a ODE represented by an ODEOperator `op` an initial condition `u0` and a time span `t0,tF`, solve
# this ODE lazily with the provided solver. I.e., return a ODESolution object (see below)
# that can be iterated to recover the pairs (u_n, t_n) at
# each computed step.
#
# The user API will be something like
#
# op = ODEOperatorMock(1.1,1.2,1.3)
#
# dt = 0.01
# nls = NLSolver()
# solver = BackwardEuler(nls,dt)
# @santiagobadia: Think about dt that can change in time, e.g., using a mutable
# array that is just dt all time steps but that could allow more complicated situations,
# e.g., time adaptivity
#
# u0 = [0,0]
# steps = solve(solve,op,u0,0,10)
#
# for (u_n, t_n) in steps
#
#   println("The solution at time $(t_n) is $(u_n)")
#
# end
#
function solve(
  solver::ODESolver,op::ODEOperator,u0::AbstractVector,t0::Real,tF::Real)
  GenericODESolution(solver,op,u0,t0,tF)
end

# I think we can create the machinery for ODE solvers easily, using more or less
# what we have above. But there is still a layer that combines it with the FE
# machinery in Gridap.
