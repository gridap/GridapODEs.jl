using Gridap
using Test
using GridapTimeStepper.ODETools
using GridapTimeStepper.TransientFETools


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
U0 = U(1.0)
ud0 = get_dirichlet_values(U0)
U1 = U(2.0)
ud1 = get_dirichlet_values(U1)
@test all(ud0 .≈ 0.5ud1)

Ut = ∂t(U)
Ut0 = Ut(0.0)
Ut1 = Ut(1.0)
utd0 = get_dirichlet_values(Ut0)
utd1 = get_dirichlet_values(Ut1)
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

t0 = 0.0
tF = 1.0
dt = 0.1

# u0 = zeros(V0)

using LineSearches: BackTracking
nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking())
nlfes = FESolver(nls)

# We create a ODE solver, e.g., Backward-Euler and solve the problem
# @santiagobadia : BackwardEuler requires a NonlinearSolver but we
# also need for FE machinery a NonlinearFESolver. What to do here?
# I don't want to create a method for every solver... How can I solve it?
be = BackwardEuler(nlfes.nls,dt)
sbt = subtypes(ODESolver)
bes = sbt[1]

# eval(Meta.parse(bes))
# (x) = 2*x

# be(1.0)
using Gridap.FESpaces: NonlinearFESolver
import GridapTimeStepper.ODETools: BackwardEuler
##
# for odes in subtypes(ODESolver)
#   (::odes)(nlfes::NonlinearFESolver,dt) = 2.0 # a(nlfes.nls,dt)
# end
# bes
#
# (::bes)(nlfes::NonlinearFESolver,dt) = 2.0 # a(nlfes.nls,dt)
# bes(nlfes,dt)
# # isa(nlfes,NonlinearFESolver)
# nbe = BackwardEuler(nlfes,dt)




nlfes



be(1.0,1.0)

sol_t = solve(be,op,u0,t0,tF)

l2(w) = w*w
h1(w) = a(w,w) + l2(w)
##

# We test it ...
for (uh_tn, tn) in sol_t
  u(x::Point) = u(x,tn)
  ∇u(x::Point) = ∇u(x,tn)

  e = u(tn) - uh_tn
  el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
  eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))
  @test el2 < tol
  @test eh1 < tol
  # writevtk(trian,"sol at time: $tn",cellfields=["u" => uh_tn])
end
