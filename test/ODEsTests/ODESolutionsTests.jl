
using Gridap
using GridapTimeStepper
using GridapTimeStepper.ODETools: BackwardEuler

include("ODEOperatorMocks.jl")

op = ODEOperatorMock(1.1,1.2,1.3)

dt = 0.01
nls = NLSolver()

solver = BackwardEuler(nls,dt)
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

solver = BackwardEuler()

struct GenericODESolution <: ODESolution
  solver::ODESolution
  op::ODEOperator
  u0::AbstractVector
  t0::Real
  tF::Real
end
