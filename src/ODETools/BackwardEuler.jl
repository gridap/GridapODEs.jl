struct BackwardEuler <: ODESolver
  nls::NonLinearSolver
  dt::Float64
end

get_step_size(self::BackwardEuler) = self.dt

function solve_step!(
  uf::AbstractVector,solver::BackwardEuler,op::ODEOperator,u0::AbstractVector,t0::Real, cache) # -> (uF,tF)

  # Build the nonlinear problem to solve at this step
  dt = solver.dt
  tf = t0+dt
  nlop = BackwardEulerNonLinearOperator(op,tf,dt,u0) # See below

  # Solve the nonlinear problem
  if (cache==nothing)
    cache = solve!(uf,solver.nls,nlop)
  else
    solve!(uf,solver.nls,nlop,cache)
  end

  # Return pair
  return (uf, tf, cache)
end

function allocate_cache(
  solver::BackwardEuler,op::ODEOperator,u0::AbstractVector,t0::Real)
  r = allocate_residual(solver.nls,x)
  J = allocate_jacobian(solver.nls,x)
  dx = copy(u0)
  (r, J, dx)
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
  vF = (x-op.u0)/op.dt
  residual!(b,op.odeop,op.tF,uF,vF)
end

# @santiagobadia : TO BE CHANGED, just a hack!!!
function fill_entries!(J,v)
  J .= zero(eltype(J))
end

function jacobian!(A::AbstractMatrix,op::BackwardEulerNonLinearOperator,x::AbstractVector)
  uF = x
  vF = (x-op.u0)/op.dt
  fill_entries!(A,zero(eltype(A)))
  jacobian_unknown!(A,op.odeop,op.tF,uF,vF)
  jacobian_unknown_t!(A,op.odeop,op.tF,uF,vF,(1/op.dt))
end

function allocate_residual(op::BackwardEulerNonLinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,x,x)
end

function allocate_jacobian(op::BackwardEulerNonLinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,x,x)
end

function zero_initial_guess(::Type{T},op::BackwardEulerNonLinearOperator) where T
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end
