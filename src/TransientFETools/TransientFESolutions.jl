"""
It represents a FE function at a set of time steps. It is a wrapper of a ODE
solution for free values combined with data for Dirichlet values. Thus, it is a
lazy iterator that computes the solution at each time step when accessing the
solution.
"""
struct TransientFESolution
  odesol::ODESolution
end

function TransientFESolution(solver::TransientFESolver,
                             op::TransientFEOperator,
                             uh0,
                             t0::Real,
                             tF::Real)
  odes = solver.odes
  ode_op = get_algebraic_operator(op)
  u0 = get_free_values(uh0)
  ode_sol = GenericODESolution(odes,ode_op,u0,t0,tF)
  TransientFESolution(ode_sol)
end

function Base.iterate(sol::TransientFESolution)
  (uf, tf), state = Base.iterate(sol.odesol)
  uf,u0,tf,op_state,cache = state
  Uh ,Uht = op_state
  uh = FEFunction(Uh,uf)
  (uh, tf), state
end

function Base.iterate(sol::TransientFESolution, state)
  uf,u0,tf,op_state,cache = state
  if tf > sol.odesol.tF
    return nothing
  end
  (uf, tf), state = Base.iterate(sol.odesol,state)
  uf,u0,tf,op_state,cache = state
  Uh ,Uht = op_state
  uh = FEFunction(Uh,uf)
  (uh, tf), state
end
