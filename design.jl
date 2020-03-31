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


# @santiagobadia Here you are defining A(t,u,u_t) in A(t,u,u_t) = 0.
abstract type ODEOperator <: GridapType end
# @santiagobadia: There is quite a lot of noise about what an ODE is, i.e.,
# whether it requires a non-singular Jacobian wrt u_t (DAE in this case) or not.
# I decided to use the term TimeStepper to circumvent this issue. In some places
# they use the term implicit ODEs to what I propose to use.

# For me the ODEOperator should be quite simple. In the Gridap FEM context,
# it will take a residual r(t,u,u_t,v), jacobian_unk(t,u,u_t,v) and
# jacobian_unk_t(t,u,u_t,v). In any case, there are some things to be
# considered when designing this object. TrialSpaces can also depend on time,
# due to time dependent boundary conditions. Also, new transient functions
# will be needed, e.g., u_D(x,t). We must think how to implement all this.

# This ones represents the value `A(t,u,v)`
function residual!(r::AbstractVector,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  @abstractmethod
end
# @santiagobadia : Agreed


function allocate_residual(op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  @abstractmethod
end
# @santiagobadia : Agreed

# This one represents `[∂A/∂u](t,u,v)`
function jacobian_unk!(j_u::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  @abstractmethod
end
# @santiagobadia : Agreed
# @santiagobadia : I think unk and unk_t or u u_t is more revealing than u v
# Not exposed to user anyway

# This one represents `[∂A/∂v](t,u,v)`
function jacobian_unk_t!(j_u_t::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  @abstractmethod
end
# @santiagobadia : Agreed

function allocate_jacobian(op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  @abstractmethod
end
# @santiagobadia : Agreed

# EXAMPLE:

# Toy linear ODE with 2 DOFs
# u_1_t = a * u_1
# u_2_t = b * u_1 + c * u_2
struct ODEOperatorMock{T<:Real} <: ODEOperator
  a::T
  b::T
  c::T
end

function residual!(r::AbstractVector,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  r .= 0
  r[1] = - u_t[1] + op.a * u[1]
  r[2] = - u_t[2] + op.b * u[1] + op.c * u[2]
  r
end

function allocate_residual(op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  zeros(2)
end

function jacobian_unk!(j_u::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  j_u .= 0
  j_u[1,1] = op.a
  j_u[2,1] = op.b
  j_u[2,2] = op.c
  j_u
end

function jacobian_unk_t!(j_u_t::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  j_u_t .= 0
  j_u_t[1,1] = -1
  j_u_t[2,2] = -1
  j_u_t
end

function allocate_jacobian(op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
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

abstract type ODESolution <: GridapType end
# @santiagobadia : Agreed
# @santiagobadia : I think this type is not a TimeStepper for Gridap, it is
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
  dt::Float64 # @santiagobadia : Previous comment about dt
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
  # @santiagobadia : Does it work so simple? Implemented for CSX matrices?
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

# I think we can create the machinery for ODE solvers easily, using more or less
# what we have above. But there is still a layer that combines it with the FE
# machinery in Gridap. 
