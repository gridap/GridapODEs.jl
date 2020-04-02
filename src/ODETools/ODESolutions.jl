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

# @santiagobadia: Getters for writing it general?
# not sure we will need more types like this one

function test_ode_solution(sol::ODESolution)
  uf = copy(sol.u0) # getter
  dt = sol.solver.dt # getter
  t0 = sol.t0
  
  current, state = Base.iterate(sol)
  uf, tf = current
  uf, u0, tf, cache = state
  @test tf==t0+dt

  current, state = Base.iterate(sol,state)
  uf, tf = current
  @test tf≈t0+2*dt
  uf, u0, tf, cache = state

  _t_n = t0
  for (u_n, t_n) in sol
    _t_n += dt
    @test t_n≈_t_n
  end

end
