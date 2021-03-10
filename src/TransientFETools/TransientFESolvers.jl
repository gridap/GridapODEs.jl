"""
A general TransientFESolver
"""
struct TransientFESolver
  odes::ODESolver
end

function solve(
  solver::TransientFESolver,op::TransientFEOperator,u0,t0::Real,tf::Real)
  TransientFESolution(solver,op,u0,t0,tf)
end

function solve(
  solver::TransientFESolver,op::TransientFEOperator,u0,v0,a0,t0::Real,tf::Real)
  TransientFESolution(solver,op,u0,v0,a0,t0,tf)
end

function test_transient_fe_solver(solver::TransientFESolver,op::TransientFEOperator,u0,t0,tf)
  solution = solve(solver,op,u0,t0,tf)
  test_transient_fe_solution(solution)
end
