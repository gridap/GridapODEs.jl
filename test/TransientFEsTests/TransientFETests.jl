using Gridap
using Test
using GridapTimeStepper.ODETools
using GridapTimeStepper.TransientFETools
using Gridap.FESpaces: get_algebraic_operator

# First, we define the transient problem
u(x,t) = (x[1] + x[2])*t
u(t::Real) = x -> u(x,t)
∇u(x,t) = VectorValue(1,1)*t
∇u(t::Real) = x -> ∇u(x,t)
import Gridap: ∇
∇(::typeof(u)) = ∇u
∇(u) === ∇u

# u(x::Point) = u(x,0.0)
# p = Point(1.0,1.0)
# u(p)
# for tn in 0:10
#   global u, ∇u
#   u(x::Point) = u(x,convert(Float64,tn))
#   ∇u(x::Point) = ∇u(x,tn)
#   @show u(p)
#   @show ∇u(p)
# end

∂tu(t) = x -> x[1]+x[2]
import GridapTimeStepper.TransientFETools: ∂t
∂t(::typeof(u)) = ∂tu
@test ∂t(u) === ∂tu

f(t) = x -> x[1]+x[2]

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

order = 1
V0 = TestFESpace(
  reffe=:Lagrangian, order=order, valuetype=Float64,
  conformity=:H1, model=model, dirichlet_tags="boundary")


U = TransientTrialFESpace(V0,u)
@test test_transient_trial_fe_space(U)

U0 = U(1.0)
ud0 = copy(get_dirichlet_values(U0))
_ud0 = get_dirichlet_values(U0)
U1 = U(2.0)
ud1 = copy(get_dirichlet_values(U1))
_ud1 = get_dirichlet_values(U1)
@test all(ud0 .≈ 0.5ud1)
all(_ud0 .≈ _ud1)

Ut = ∂t(U)
Ut0 = Ut(0.0)

using Gridap.FESpaces: TrialFESpace!
TrialFESpace!(Ut0,u(0))

Ut1 = Ut(1.0)
utd0 = copy(get_dirichlet_values(Ut0))
utd1 = copy(get_dirichlet_values(Ut1))
@test all(utd0 .== utd1)
@test all(utd1 .== ud0)

trian = Triangulation(model)
degree = 2
quad = CellQuadrature(trian,degree)

a(u,v) = ∇(v)*∇(u)
b(v,t) = v*f(t)

# Next, we create the transient and steady terms
res(t,u,ut,v) = a(u,v) + ut*v - b(v,t)
jac(t,u,ut,du,v) = a(du,v)
jac_t(t,u,ut,dut,v) = dut*v

# TransientFETerm or FETerm, what do we prefer?
t_Ω = FETerm(res,jac,jac_t,trian,quad)
# We create the transient operator
op = TransientFEOperator(U,V0,t_Ω)
U0 = U(0.0)
u0 = interpolate_everywhere(U0,0.0)
_u0 = get_free_values(u0)
@test test_transient_fe_operator(op,u0)

u0 = u(0.0)
t0 = 0.0
tF = 1.0
dt = 0.1

# u0 = zeros(V0)
ls = LUSolver()
# using LineSearches: BackTracking
tol = 1.0
maxiters = 20
using Gridap.Algebra: NewtonRaphsonSolver
# nls = NewtonRaphsonSolver(ls,tol,maxiters)
nls = NLSolver(ls;show_trace=true,method=:newton) #linesearch=BackTracking())
odes = BackwardEuler(nls,dt)
solver = TransientFESolver(odes) # Return a specialization of TransientFESolver
@test test_transient_fe_solver(solver,op,u0,t0,tF)

odeop = get_algebraic_operator(op)
sol_ode_t = solve(odes,odeop,_u0,t0,tF)
# test_ode_solution(sol_ode_t)
# state = iterate(sol_ode_t)
# sol = sol_ode_t
#
# uf = copy(sol.u0)
# u0 = copy(sol.u0)
# cache = nothing
# t0 = sol.t0
# op_state = allocate_state(sol.op)
#
# # Solve step
# # uf, tf, op_state, cache = solve_step!(uf,sol.solver,sol.op,u0,t0,op_state,cache)
#
# # Update
# u0 .= uf
# state = (uf,u0,tf,cache)
#
#
#
# sol = sol_ode_t
#
# _t_n = t0
# Base.iterate(sol_ode_t)
# for (u_n, t_n) in sol_ode_t
#   global _t_n
#   _t_n += dt
#   @test t_n≈_t_n
# end

# @santiagobadia : Now it is time to check the FE problem !!!

# ##
# # Base.iterate(sol_ode_t)
# uf = copy(sol.u0)
# u0 = copy(sol.u0)
# cache = nothing
# t0 = sol.t0
# op_state = allocate_state(sol.op)
#
# solver = odes
#
#
#
# op = sol.op
# # uf, tf, op_state, cache = solve_step!(uf,sol.solver,sol.op,u0,t0,op_state,cache)
# # function solve_step!(
# #   uf::AbstractVector,solver::BackwardEuler,op::ODEOperator,u0::AbstractVector,t0::Real,op_state,cache) # -> (uF,tF)
# dt = solver.dt
# tf = t0+dt
# update_state!(op_state,op,tf)
# using GridapTimeStepper.ODETools: BackwardEulerNonlinearOperator
# nlop = BackwardEulerNonlinearOperator(op,tf,dt,u0,op_state) # See below
# # cache =
# solve!(uf,solver.nls,nlop)
# using GridapTimeStepper.ODETools: allocate_residual
# allocate_residual(nlop,u0)
# # @show nlop.op_state
# allocate_residual(nlop.odeop,u0)
#
# # Solve the nonlinear problem
# # if (cache==nothing)
# # else
# #   solve!(uf,solver.nls,nlop,cache)
# # end
#
#
#
# # uf, tf, op_state, cache =
# solve_step!(uf,sol.solver,sol.op,u0,t0,op_state,cache)
# # @test test_ode_operator(odeop,0.0,_u0,_u0)
#
# # test_ode_solution(sol_ode_t)
#
#
# # @test test_ode_solver(odes,odeop,_u0,t0,tF)
#
#
# #
#









# sol_t = solve(solver,op,u0,t0,tF)

l2(w) = w*w
h1(w) = a(w,w) + l2(w)
##

# We test it ...
# for (uh_tn, tn) in sol_t
#   u(x::Point) = u(x,tn)
#   ∇u(x::Point) = ∇u(x,tn)
#
#   e = u(tn) - uh_tn
#   el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
#   eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))
#   @test el2 < tol
#   @test eh1 < tol
#   # writevtk(trian,"sol at time: $tn",cellfields=["u" => uh_tn])
# end
