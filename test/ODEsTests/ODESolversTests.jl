# module ODESolversTests

using GridapTimeStepper
using GridapTimeStepper.ODETools: GenericODESolution
using GridapTimeStepper.ODETools: BackwardEuler
using GridapTimeStepper.ODETools: BackwardEulerNonLinearOperator
using GridapTimeStepper.ODETools: solve!
using GridapTimeStepper
using GridapTimeStepper.ODETools
using Gridap
using Test

# using Gridap.Algebra: residual, jacobian

include("ODEOperatorMocks.jl")

op = ODEOperatorMock(1.0,0.0,1.0)

include("ODESolverMocks.jl")

t0 = 0.0
tf = 1.0
dt = 0.1
u0 = ones(2)*2

# NonlinearOperator tests

sop = OperatorMock(op,tf,dt,u0)

x = zero_initial_guess(eltype(u0),sop)
x .+= 1.0
r = allocate_residual(sop,x)
J = allocate_jacobian(sop,x)
residual!(r,sop,x)
jacobian!(J,sop,x)
@test all(r .== [ -11.0 -11.0])
@test all(J .== [ 9.0 0.0; 0.0 9.0])
_r = residual(sop,x)
_J = jacobian(sop,x)
@test all(_r .== [ -11.0 -11.0])
@test all(_J .== [ 9.0 0.0; 0.0 9.0])

# NLSolver tests

nls = NLSolverMock()
cache = solve!(x,nls,sop)
r, J, dx = cache
@test all(r.==_r)
@test all(J.==_J)
@test all(dx.≈11/9)
@test all(x.≈1+11/9)

#ODESolver tests

odesol = ODESolverMock(nls,dt)
uf = copy(u0)
uf.=1.0
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,nothing)
uf
@test tf==t0+dt
@test all(uf.≈x)

# ODESolutions

tF = 10.0
sol = GenericODESolution(odesol,op,u0,t0,tF)
current, state = Base.iterate(sol)
uf, tf = current
@test tf==t0+dt
@test all(uf.≈x)

# BackwardEulerNonLinearOperator tests

tf = t0+dt
sop = BackwardEulerNonLinearOperator(op,tf,dt,u0) # See below
x = zero_initial_guess(eltype(u0),sop)
x .+= 1.0
r = allocate_residual(sop,x)
J = allocate_jacobian(sop,x)
residual!(r,sop,x)
jacobian!(J,sop,x)
@test all(r .== [ -11.0 -11.0])
@test all(J .== [ 9.0 0.0; 0.0 9.0])
_r = residual(sop,x)
_J = jacobian(sop,x)
@test all(_r .== [ -11.0 -11.0])
@test all(_J .== [ 9.0 0.0; 0.0 9.0])


# BackwardEuler tests

# @santiagobadia : I have a question, how to use an existing NL Solver
# One from Julia package? One from Gridap? Dif's? @fverdugo help on this
odesol = BackwardEuler(nls,dt)
uf = copy(u0)
uf.=1.0
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,nothing)
uf
@test tf==t0+dt
@test all(uf.≈1+11/9)

@test test_ode_solver(odesol,op,u0,t0,uf,tf)

# end #module
