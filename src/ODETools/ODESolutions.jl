# Represents a lazy iterator over all solution in a time interval
abstract type ODESolution <: GridapType end

# First time step
function Base.iterate(u::ODESolution) # (u0,t0)-> (uf,tf) or nothing
  @abstractmethod
end

# Following time steps
function Base.iterate(u::ODESolution,state) # (u0,t0)-> (uf,tf) or nothing
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

  # Solve step
  uf, tf, cache = solve_step!(uf,sol.solver,sol.op,u0,t0,cache)

  # Update
  u0 .= uf
  state = (uf,u0,tf,cache)

  return (uf, tf), state
end

function Base.iterate(sol::GenericODESolution, state)

  uf,u0,t0,cache = state

  if t0 > sol.tF
    return nothing
  end

  # Solve step
  uf, tf, cache = solve_step!(uf,sol.solver,sol.op,u0,t0,cache)

  # Update
  u0 .= uf
  state = (uf,u0,tf,cache)

  return (uf, tf), state
end
