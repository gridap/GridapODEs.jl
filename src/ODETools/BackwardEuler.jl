struct BackwardEuler <: ODESolver
  nls::NonLinearSolver
  dt::Float64
end

function solve_step!(
  uF::AbstractVector,solver::BackwardEuler,op::ODEOperator,u0::AbstractVector,t0::Real, cache) # -> (uF,tF)

  # Build the non-linear problem to solve at this step
  dt = solver.dt
  tF = t0+dt
  nlop = BackwardEulerNonLinearOperator(op,tF,dt,u0) # See below

  # Solve the nonlinear problem
  # uF, cache = solve!(uF,solver.nls,nlop,cache) # TODO reuse the cache
  uF, cache = (u0*tF,cache)

  # Return pair
  return (uF, tF, cache)
end

function allocate_cache(
  solver::BackwardEuler,op::ODEOperator,u0::AbstractVector,t0::Real)
  nothing
end

# Struct representing the nonlinear algebraic problem to be solved at a given step
struct BackwardEulerNonLinearOperator <: NonLinearOperator
  odeop::ODEOperator
  tF::Float64
  dt::Float64
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
  fill_entries!(A,zero(eltype(A)))
  jacobian_unknown!(A,op.odeop,op.tF,uF,vF)
  jacobian_unknown_t!(A,op.odeop,op.tF,uF,vF,(1/op.dt))
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
