struct BackwardEuler <: ODESolver
  nls::NonlinearSolver
  dt::Float64
end

function solve_step!(
  uf::AbstractVector,solver::BackwardEuler,op::ODEOperator,u0::AbstractVector,t0::Real,op_state,cache) # -> (uF,tF)

  # Build the nonlinear problem to solve at this step
  dt = solver.dt
  tf = t0+dt
  update_state!(op_state,op,tf)
  nlop = BackwardEulerNonlinearOperator(op,tf,dt,u0,op_state) # See below


  # Solve the nonlinear problem
  if (cache==nothing)
    cache = solve!(uf,solver.nls,nlop)
  else
    solve!(uf,solver.nls,nlop,cache)
  end

  # Return pair
  return (uf,tf,op_state,cache)
end

# Struct representing the nonlinear algebraic problem to be solved at a given step
struct BackwardEulerNonlinearOperator <: NonlinearOperator
  odeop::ODEOperator
  tF::Float64
  dt::Float64
  u0::AbstractVector
  op_state
end

function residual!(b::AbstractVector,op::BackwardEulerNonlinearOperator,x::AbstractVector)
  uF = x
  vF = (x-op.u0)/op.dt
  residual!(b,op.odeop,op.tF,uF,vF,op.op_state)
end

function jacobian!(A::AbstractMatrix,op::BackwardEulerNonlinearOperator,x::AbstractVector)
  uF = x
  vF = (x-op.u0)/op.dt
  z = zero(eltype(A))
  fill_entries!(A,z)
  jacobian!(A,op.odeop,op.tF,uF,vF,op.op_state)
  jacobian_t!(A,op.odeop,op.tF,uF,vF,(1/op.dt),op.op_state)
end

function allocate_residual(op::BackwardEulerNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,x,op.op_state)
end

function allocate_jacobian(op::BackwardEulerNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,x,op.op_state)
end

function zero_initial_guess(::Type{T},op::BackwardEulerNonlinearOperator) where T
  x0 = similar(op.u0,T)
  fill!(x0,zero(T))
  x0
end
