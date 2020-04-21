# module StokesTaylorHoodTests

using Test
using Gridap
import Gridap: ∇
using GridapODEs.ODETools
using GridapODEs.TransientFETools
import GridapODEs.TransientFETools: ∂t

using LinearAlgebra: tr

u(x) = VectorValue( x[1]^2 + 2*x[2]^2, -x[1]^2 )
∇u(x) = TensorValue( 2*x[1], 4*x[2], -2*x[1], zero(x[1]) )
Δu(x) = VectorValue( 6, -2 )

p(x) = x[1]-x[2]
∇p(x) = VectorValue(1,-1)

f(x) = -Δu(x) + ∇p(x)
g(x) = tr(∇u(x))

∇(::typeof(u)) = ∇u
∇(::typeof(p)) = ∇p

u(t::Real) = x -> u(x)

domain = (0,2,0,2)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",[1,2,5])
add_tag_from_tags!(labels,"neumann",[6,7,8])

order = 2

V = TestFESpace(
model=model,
order=order,
reffe=:Lagrangian,
labels=labels,
valuetype=VectorValue{2,Float64},
dirichlet_tags="boundary",
# dirichlet_tags="dirichlet",
# dof_space=ref_st,
conformity=:H1)

Q = TestFESpace(
model=model,
order=order-1,
reffe=:Lagrangian,
valuetype=Float64,
# dof_space=ref_st,
conformity=:H1,
constraint=:zeromean)

U = TrialFESpace(V,u)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

trian = get_triangulation(model)
degree = order*2
quad = CellQuadrature(trian,degree)

btrian = BoundaryTriangulation(model,labels,"neumann")
bdegree = order
bquad = CellQuadrature(btrian,bdegree)
n = get_normal_vector(btrian)

function a(x,y)
  u,p = x
  v,q = y
  inner(∇(v),∇(u)) - (∇*v)*p + q*(∇*u)
end

function l(y)
  v,q = y
  v*f + q*g
end

function l_Γb(y)
  v,q = y
  v*(n*∇u) - (n*v)*p
end

t_Ω = AffineFETerm(a,l,trian,quad)
# t_Γb = FESource(l_Γb,btrian,bquad)

# op = AffineFEOperator(X,Y,t_Ω,t_Γb)
op = AffineFEOperator(X,Y,t_Ω)

uh, ph = solve(op)

eu = u - uh
ep = p - ph

l2(v) = v*v
h1(v) = v*v + inner(∇(v),∇(v))

eu_l2 = sqrt(sum(integrate(l2(eu),trian,quad)))
eu_h1 = sqrt(sum(integrate(h1(eu),trian,quad)))
ep_l2 = sqrt(sum(integrate(l2(ep),trian,quad)))

tol = 1.0e-9
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

θ = 0.5


u(x,t) = u(x)
∂tu(t) = 0.0
∂tu(x,t) = ∂tu(t)(x)
u(t::Real) = x -> u(x,t)

∂t(::typeof(u)) = ∂tu
∇(::typeof(u)) = ∇u

U = TransientTrialFESpace(V,u)
Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

function res(t,x,xt,y)
  u,p = x
  ut,pt = xt
  v,q = y
  a(x,y) + 0.0*ut*v - l(y)
end

function jac(t,x,xt,dx,y)
  a(dx,y)
end

function jac_t(t,x,xt,dxt,y)
  dut,dpt=dxt
  v,q = y
  dut*v
end

t_Ω = FETerm(res,jac,jac_t,trian,quad)
t_Γb = FESource(l_Γb,btrian,bquad)
# op = TransientFEOperator(X,Y,t_Ω,t_Γb)
op = TransientFEOperator(X,Y,t_Ω)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(U0,u(0.0))
ph0 = interpolate_everywhere(P,0.0)

X0 = X(0.0)
xh0 = Gridap.MultiField.MultiFieldFEFunction(X0,[uh0,ph0])

ls = LUSolver()
# using Gridap.Algebra: NewtonRaphsonSolver
# nls = NLSolver(ls;show_trace=true,method=:newton) #linesearch=BackTracking())
odes = ThetaMethod(ls,dt,θ)
solver = TransientFESolver(odes) # Return a specialization of TransientFESolver

sol_t = solve(solver,op,xh0,t0,tF)

l2(w) = w*w


tol = 1.0e-6
_t_n = t0

result = Base.iterate(sol_t)


for (xh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  uh_tn = xh_tn.blocks[1]
  e = u(tn) - uh_tn
  el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
  @test el2 < tol
end

# end #module
