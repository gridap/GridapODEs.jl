# module ODESolversTests

using GridapODEs
using GridapODEs.ODETools: GenericODESolution
using GridapODEs.ODETools: BackwardEuler
using GridapODEs.ODETools: ThetaMethodNonlinearOperator
using GridapODEs.ODETools: solve!
using GridapODEs
using GridapODEs.ODETools
using Gridap
using Test

# using Gridap.Algebra: residual, jacobian

include("ODEOperatorMocks.jl")

op = ODEOperatorMock{Float64,Constant}(1.0,0.0,1.0)

include("ODESolverMocks.jl")

t0 = 0.0
tf = 1.0
dt = 0.1
u0 = ones(2)*2

# NonlinearOperator tests

sop = OperatorMock(op,tf,dt,u0)
isa(sop,NonlinearOperator)

ode_cache = allocate_cache(op)

x = zero_initial_guess(sop)
x .+= 1.0
isa(sop,OperatorMock)
isa(x,AbstractVector)
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

# BackwardEulerNonlinearOperator tests

tf = t0+dt
vf = copy(u0)
sop = ThetaMethodNonlinearOperator(op,tf,dt,u0,ode_cache,vf) # See below
x = zero_initial_guess(sop)
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

ls = LUSolver()
odesol = BackwardEuler(ls,dt)
uf = copy(u0)
uf.=1.0
cache = nothing
# Juno.@enter solve_step!(uf,odesol,op,u0,t0,cache)
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all(uf.≈1+11/9)

@test test_ode_solver(odesol,op,u0,t0,tf)
test_ode_solver(odesol,op,u0,t0,tf)

# Affine and nonlinear solvers

op = ODEOperatorMock{Float64,Nonlinear}(1.0,0.0,1.0)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all(uf.≈1+11/9)

op = ODEOperatorMock{Float64,Affine}(1.0,0.0,1.0)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all(uf.≈1+11/9)


# end #module
