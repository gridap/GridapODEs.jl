using Gridap
using Test

# First, we define the transient problem
u(t) = x -> (x[1] + x[2])*t
∇u(t) = x -> VectorValue(1,1)*t
# @santiagobadia : Better way to do this?
# u(t,x) don't think it will be possible
import Gridap: ∇
∇(::typeof(u)) = ∇u
∇(u) === ∇u
# @santiagobadia: It is not going to work internally...
# It is not the gradient of u but its result for a given t what we need
# to link to its gradient... using return_type...

∂tu = x[1]+x[2]
import GridapTimeStepper: ∂t
∂t(::typeof(u)) = ∂tu
∂t(u) === ∂tu

f(t) = x -> x[1]+x[2]

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

order = 1
V0 = TestFESpace(
  reffe=:Lagrangian, order=order, valuetype=Float64,
  conformity=:H1, model=model, dirichlet_tags="boundary")

U = TransientTrialFESpace(V0,t->u(t))

trian = Triangulation(model)
degree = 2
quad = CellQuadrature(trian,degree)

a(u,v) = ∇(v)*∇(u)
b(v) = v*f

t_Ω = AffineFETerm(a,b,trian,quad)
op = AffineFEOperator(U,V0,t_Ω)

uh = solve(op)

# Next, we create the transient and steady terms

res(u,ut,v) = a(u,v) + ut*v - b(v)
jac(u,ut,du,v) = a(du,v)
jac_t(u,ut,dut,v) = dut*v

t_Ω = FETerm(res,jac,jac_t,trian,quad)
# We create the transient operator
op = TransientFEOperator(V,t -> U,t_Ω)
# I don't think we need here a space for the time derivative, using ∂t
# @santiagobadia : We will need to define a ∂t method a la ∇ to compute
# internally the time derivative of a Dirichlet boundary
# condition

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
  e = u(tn) - uh_tn
  el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
  eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))
  @test el2 < tol
  @test eh1 < tol
  # writevtk(trian,"sol at time: $tn",cellfields=["u" => uh_tn])
end
