# Like we have done in the spatial FE discretization, I would distinguish between the
# algebraic and FE counterparts and between the problem to be solved and the method to solve the problem.
#
# Disclaimer: The names of functions and types I use here are just an initial guess. They can be improved a lot.
#
# For the moment, we can focus in the algebraic part.
# From this, designing the FE counterpart will be easy.
#
# First ingredient: object describing the algebraic problem to solve.
# We can follow the notation A(t,u,u_t)=0 you have proposed to represent the ODE.
# For simplicity, I use here `A(t,u,v)=0`.
# The first we need is an abstract type representing the operator A:
#

abstract type ODEOperator <: GridapType end

# This ones represents the value `A(t,u,v)`
function residual!(r::AbstractVector,op::ODEOperator,t::Real,u::AbstractVector,v::AbstractVector)
  @abstractmethod
end

function allocate_residual(op::ODEOperator,t::Real,u::AbstractVector,v::AbstractVector)
  @abstractmethod
end

# This one represents `[∂A/∂u](t,u,v)`
# Better function name?
function jacobian_u!(j_u::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,v::AbstractVector)
  @abstractmethod
end

# This one represents `[∂A/∂v](t,u,v)`
# Better function name?
function jacobian_v!(j_v::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,v::AbstractVector)
  @abstractmethod
end

function allocate_jacobian(op::ODEOperator,t::Real,u::AbstractVector,v::AbstractVector)
  @abstractmethod
end

# EXAMPLE:

# Toy linear ODE with 2 DOFs
# v_1 = a * u_1
# v_2 = b * u_1 + c * u_2 
struct ODEOperatorMock{T<:Real} <: ODEOperator
  a::T
  b::T
  c::T
end

function residual!(r::AbstractVector,op::ODEOperatorMock,t::Real,u::AbstractVector,v::AbstractVector)
  r .= 0
  r[1] = - v[1] + op.a * u[1]
  r[2] = - v[2] + op.b * u[1] + op.c * u[2]
  r
end

function allocate_residual(op::ODEOperatorMock,t::Real,u::AbstractVector,v::AbstractVector)
  zeros(2)
end

function jacobian_u!(j_u::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,v::AbstractVector)
  j_u .= 0
  j_u[1,1] = op.a
  j_u[2,1] = op.b
  j_u[2,2] = op.c
  j_u
end

function jacobian_v!(j_v::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,v::AbstractVector)
  j_v .= 0
  j_v[1,1] = -1
  j_v[2,2] = -1
  j_v
end

function allocate_jacobian(op::ODEOperatorMock,t::Real,u::AbstractVector,v::AbstractVector)
  zeros(2,2)
end

# Now, we need an abstract type representing a numerical discretization scheme for the ODE
#
abstract type ODESolver <: GridapType end

# The solver is defined by how a single step is solved
function solve_step!(
  uF::AbstractVector,solver::ODESolver,op::ODEOperator,u0::AbstractVector,t0::Real) # -> (uF,tF)
  @abstractmethod
end


# Assuming that solve_step! is defined, we can define the following solve function for all ODESolver objects:
# Given an ODE represented by an ODEOperator `op` an initial condition `u0` and a time span `t0,tF`, solve 
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

## This is what you have called TimeStepper
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

function Base.iterate(sol::GenericODESolution)

  uF = copy(sol.u0)
  u0 = copy(sol.u0)

  # Solve step
  uF, tF = solve_step!(uF,sol.op,u0,sol.t0)

  # Update
  u0 .= uF
  state = (uF,u0,tF)

  return (uf, tF), state
end

function Base.iterate(sol::GenericODESolution, state)

  uF,u0,t0 = state

  if t0 > op.tF
    return nothing
  end

  # Solve step
  uF, tF = solve_step!(uF,sol.op,u0,t0)

  # Update
  u0 .= uF
  state = (uF,u0,tF)

  return (uf, tF), state
end

# EXAMPLE
#  Vanilla backward Euler

struct BackwardEuler <: ODESolver
  nls::NonLinearSolver
  dt::Float64
end

function solve_step!(
  uF::AbstractVector,solver::BackwardEuler,op::ODEOperator,u0::AbstractVector,t0::Real) # -> (uF,tF)

  # Build the non-linear problem to solve at this step
  dt = op.dt
  tF = t0+dt
  nlop = BackwardEulerNonLinearOperator(op,tF,dt,u0) # See below

  # Solve the nonlinear problem
  uF, cache = solve!(uF,solver.nls,nlop) # TODO reuse the cache

  # Return pair
  return (uF, tF)
end

# Struct representing the non-linear algebraic problem to be solved at a given step
struct BackwardEulerNonLinearOperator <: NonLinearOperator
  odeop::ODEOperator
  tF::Float64
  dt = Float64
  u0::AbstractVector
end

function residual!(b::AbstractVector,op::BackwardEulerNonLinearOperator,x::AbstractVector)
  uF = x
  vF = (x-u0)/op.dt
  residual!(b,op.odeop,op.tF,uF,vF)
end

function jacobian!(A::AbstractMatrix,op::BackwardEulerNonLinearOperator,x::AbstractVector)
  uF = x
  vF = (x-u0)/op.dt
  A_u = copy(A) # TODO avoid this
  A_v = copy(A) # TODO avoid this
  jacobian_u!(A_u,op.odeop,op.tF,uF,vF)
  jacobian_v!(A_v,op.odeop,op.tF,uF,vF)
  A .= A_u + (1/op.dt)*A_v
end

function allocate_residual(op::BackwardEulerNonLinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,op.odeop.tF,x,x)
end

function allocate_jacobian(op::BackwardEulerNonLinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,op.odeop.tF,x,x)
end

function zero_initial_guess(::Type{T},op::BackwardEulerNonLinearOperator) where T
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end


