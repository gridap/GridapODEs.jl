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

u(x::Point) = u(x,0.0)
p = Point(1.0,1.0)
u(p)
for tn in 0:10
  global u, ∇u
  u(x::Point) = u(x,convert(Float64,tn))
  ∇u(x::Point) = ∇u(x,tn)
  @show u(p)
  @show ∇u(p)
end

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

##
U = TransientTrialFESpace(V0,u)
U0 = U(0.0)
get_dirichlet_values(U0)
Ut = ∂t(U)
Ut0 = Ut(0.0)
get_dirichlet_values(Ut0)
i##
trian = Triangulation(model)
degree = 2
quad = CellQuadrature(trian,degree)

a(u,v) = ∇(v)*∇(u)
b(v,t) = v*f(t)

# Next, we create the transient and steady terms

res(t,u,ut,v) = a(u,v) + ut*v - b(v,t)
jac(t,u,ut,du,v) = a(du,v)
jac_t(t,u,ut,dut,v) = dut*v

t_Ω = TransientFETerm(res,jac,jac_t,trian,quad)

# We create the transient operator
op = TransientFEOperator(V,U,t_Ω)

t0 = 0.0
tF = 1.0
dt = 0.1

u0 = zeros(V)

using LineSearches: BackTracking
nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking())
nlfes = FESolver(nls)

# We create a ODE solver, e.g., Backward-Euler and solve the problem
be = BackwardEuler(nlfes,dt)
sol_t = solve(be,op,u0,t0,tF)

l2(w) = w*w
h1(w) = a(w,w) + l2(w)

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
