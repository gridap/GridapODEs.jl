
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

# tester

function test_ode_solution(sol::ODESolution)
  for (u_n,t_n) in sol
    @test isa(t_n,Real)
    @test isa(u_n,AbstractVector)
  end
  true
end

# Specialization

struct GenericODESolution <: ODESolution
  solver::ODESolver
  op::ODEOperator
  u0::AbstractVector
  t0::Real
  tF::Real
end

function Base.iterate(sol::GenericODESolution)

  uf = copy(sol.u0)
  u0 = copy(sol.u0)
  t0 = sol.t0

  # Solve step
  uf, tf, cache = solve_step!(uf,sol.solver,sol.op,u0,t0)

  # Update
  u0 .= uf
  state = (uf,u0,tf,cache)

  return (uf, tf), state
end

function Base.iterate(sol::GenericODESolution, state)

  uf,u0,t0,cache = state

  if t0 >= sol.tF - ϵ
    return nothing
  end

  # Solve step
  uf, tf, cache = solve_step!(uf,sol.solver,sol.op,u0,t0,cache)

  # Update
  u0 .= uf
  state = (uf,u0,tf,cache)

  return (uf, tf), state
end

# Specialization for a 2nd order ODE

struct Generic2ndOrderODESolution <: ODESolution
  solver::ODESolver
  op::ODEOperator
  u0::AbstractVector
  v0::AbstractVector
  a0::AbstractVector
  t0::Real
  tF::Real
end

function Base.iterate(sol::Generic2ndOrderODESolution)

  uf = copy(sol.u0)
  u0 = copy(sol.u0)
  v0 = copy(sol.v0)
  a0 = copy(sol.a0)
  t0 = sol.t0

  # Solve step
  uf, tf, cache = solve_step!(uf,sol.solver,sol.op,u0,v0,a0,t0)

  # Update
  vf, af, = cache
  u0 .= uf
  v0 .= vf
  a0 .= af
  state = (uf,u0,v0,a0,tf,cache)

  return (uf, tf), state
end

function Base.iterate(sol::Generic2ndOrderODESolution, state)

  uf,u0,v0,a0,t0,cache = state

  if t0 >= sol.tF - ϵ
    return nothing
  end

  # Solve step
  uf, tf, cache = solve_step!(uf,sol.solver,sol.op,u0,v0,a0,t0,cache)

  # Update
  vf, af, = cache
  u0 .= uf
  v0 .= vf
  a0 .= af
  state = (uf,u0,v0,a0,tf,cache)

  return (uf, tf), state
end
