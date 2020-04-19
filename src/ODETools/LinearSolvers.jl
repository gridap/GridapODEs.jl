struct LinearSolverCache <: GridapType
  A::AbstractMatrix
  b::AbstractVector
  ns#::NumericalSetup
end

function solve!(x::AbstractVector,nls::NonlinearSolver,op::NonlinearOperator,cache::Nothing)
  solve!(x,nls,op)
end

function solve!(x::AbstractVector,
                ls::LinearSolver,
                op::NonlinearOperator,
                cache::Nothing)
  x .= zero(eltype(x))
  b = residual(op, x)
  A = jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  broadcast!(*,b,b,-1)
  solve!(x,ns,b)
  LinearSolverCache(A,b,ns)
end

function solve!(x::AbstractVector,
                ls::LinearSolver,
                op::NonlinearOperator,
                cache)
  x .= zero(eltype(x))
  b = cache.b
  A = cache.A
  ns = cache.ns
  # if !isbct
  # we must recalculate it
  residual!(b, op, x)
  # end
  # if !isAct
  numerical_setup!(ns,A)
  # end
  broadcast!(*,b,b,-1)
  solve!(x,ns,b)
  cache
end

# function solve!(x::AbstractVector,
#                 ls::LinearSolver,
#                 op::NonlinearOperator,
#                 cache,
#                 isbct::Bool,
#                 isAct::Bool)
#   x .= zero(eltype(x))
#   b = cache.b
#   A = cache.A
#   ns = cache.ns
#   # if !isbct
#   # we must recalculate it
#   residual!(b, op, x)
#   # end
#   if !isAct
#     numerical_setup!(ns,A)
#   end
#   broadcast!(*,b,b,-1)
#   solve!(x,ns,b)
#   cache
# end
