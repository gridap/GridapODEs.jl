abstract type ODESolution <: GridapType end
# @santiagobadia : Agreed
# @santiagobadia : I think this type is not the full TimeStepper for Gridap, it is
# at a lower level, since it does not work with the FE machinery. The
# TimeStepper should provide FEFunction

# First time step
function Base.iterate(u::ODESolution) # -> (t_n,u_n) or nothing
  @abstractmethod
end

# Following time steps
function Base.iterate(u::ODESolution,state) # -> (t_n,u_n) or nothing
  @abstractmethod
end

# Do we want it to be type stable?
struct GenericODESolution <: ODESolution
  solver::ODESolution
  op::ODEOperator
  u0::AbstractVector
  t0::Real
  tF::Real
end

function Base.iterate(sol::GenericODESolution)

  uF = copy(sol.u0)
  u0 = copy(sol.u0)

  # Solve step
  uF, tF, cache = solve_step!(uF,sol.op,u0,sol.t0)

  # Update
  u0 .= uF
  state = (uF,u0,tF,cache)

  return (uf, tF), state
end

function Base.iterate(sol::GenericODESolution, state)

  uF,u0,t0,cache = state

  if t0 > op.tF
    return nothing
  end

  # Solve step
  uF, tF = solve_step!(uF,sol.op,u0,t0,cache)

  # Update
  u0 .= uF
  state = (uF,u0,tF,cache)

  return (uf, tF), state
end
