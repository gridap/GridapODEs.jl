"""
A general TransientFESolver
"""
struct TransientFESolver
  odes::ODESolver
end

function solve_step!(
  uF,solver::TransientFESolver,op::TransientFEOperator,u0,t0::Real,op_state,cache) # -> (uF,tF,cache)
  @abstractmethod
end

function solve(
  solver::TransientFESolver,op::TransientFEOperator,u0,t0::Real,tf::Real)
  TransientFESolution(solver,op,u0,t0,tf)
end

function test_transient_fe_solver(solver::TransientFESolver,op::TransientFEOperator,u0,t0,tf)
  solution = solve(solver,op,u0,t0,tf)
  # @santiagobadia : Do it!
  # test_transient_fe_solution(solution)
end
