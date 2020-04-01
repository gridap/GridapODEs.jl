abstract type ODESolution <: GridapType end

# First time step
function Base.iterate(u::ODESolution) # -> (t_n,u_n) or nothing
  @abstractmethod
end

# Following time steps
function Base.iterate(u::ODESolution,state) # -> (t_n,u_n) or nothing
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

  uF = copy(sol.u0)
  u0 = copy(sol.u0)

  # Solve step
  uF, tF, cache = solve_step!(uF,sol.solver,sol.op,u0,sol.t0)

  # Update
  u0 .= uF
  state = (uF,u0,tF,cache)

  return (uF, tF), state
end

function Base.iterate(sol::GenericODESolution, state)

  uF,u0,t0,cache = state

  if t0 > sol.tF
    return nothing
  end

  # Solve step
  uF, tF = solve_step!(uF,sol.solver,sol.op,u0,t0,cache)

  # Update
  u0 .= uF
  state = (uF,u0,tF,cache)

  return (uF, tF), state
end
