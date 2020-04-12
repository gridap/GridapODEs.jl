# Represents a lazy iterator over all solution in a time interval
"""
It represents the solution of a ODE at a given time interval. It is a lazy implementation,
i.e., the object is an iterator that computes the solution at each time step
when accessing the solution at each time step.
"""
abstract type ODESolution <: GridapType end

# First time step
function iterate(u::ODESolution) # (u0,t0)-> (uf,tf) or nothing
  @abstractmethod
end

# Following time steps
function iterate(u::ODESolution,state) # (u0,t0)-> (uf,tf) or nothing
  @abstractmethod
end

struct GenericODESolution <: ODESolution
  solver::ODESolver
  op::ODEOperator
  u0::AbstractVector
  t0::Real
  tF::Real
end
# @santiagobadia: Do we want it to be type stable?

function Base.iterate(sol::GenericODESolution)

  uf = copy(sol.u0)
  u0 = copy(sol.u0)
  cache = nothing
  t0 = sol.t0
  op_cache = allocate_cache(sol.op)

  # Solve step
  uf, tf, op_cache, cache = solve_step!(uf,sol.solver,sol.op,u0,t0,op_cache,cache)

  # Update
  u0 .= uf
  state = (uf,u0,tf,op_cache,cache)

  return (uf, tf), state
end

function Base.iterate(sol::GenericODESolution, state)

  uf,u0,t0,op_cache,cache = state

  if t0 > sol.tF
    return nothing
  end

  # Solve step
  uf, tf, op_cache, cache = solve_step!(uf,sol.solver,sol.op,u0,t0,op_cache,cache)

  # Update
  u0 .= uf
  state = (uf,u0,tf,op_cache,cache)

  return (uf, tf), state
end

function test_ode_solution(sol::GenericODESolution)
  for (u_n,t_n) in sol
    @test isa(t_n,Real)
    @test isa(u_n,AbstractVector)
  end
  true
end
