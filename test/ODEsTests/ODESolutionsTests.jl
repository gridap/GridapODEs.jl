
using Gridap
using GridapTimeStepper
using GridapTimeStepper.ODETools

include("ODEOperatorMocks.jl")

include("ODESolverMocks.jl")

op = ODEOperatorMock(1.1,1.2,1.3)

dt = 0.01

nls = NLSolverMock()

solver = ODESolverMock(nls,dt)
u0 = [0,0]

t0 = 0.0
tF = 10.0
steps = solve(solver,op,u0,t0,tF)
sol = steps

# uF = copy(u0)
# cache = allocate_cache(solver,op,u0,t0)
# nlop = OperatorMock(op,tF,dt,u0)
# uF, cache = solve!(uF,solver.nls,nlop,cache) # TODO reuse the cache
# uF, tF = solve_step!(uF,sol.solver,sol.op,u0,t0,cache)
#
for (u_n, t_n) in steps

  println("The solution at time $(t_n) is $(u_n)")

end
