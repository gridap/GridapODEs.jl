# module QuadraticTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator

# using GridapODEs.ODETools: ThetaMethodLinear

import Gridap: ∇
import GridapODEs.TransientFETools: ∂t

θ = 1.0

u(x,t) = VectorValue(x[1],x[2])*t
u(t::Real) = x -> u(x,t)
∂tu(t) = x -> VectorValue(x[1],x[2])
∂tu(x,t) = ∂tu(t)(x)
∂t(::typeof(u)) = ∂tu

f(t) = x -> ∂t(u)(t)(x)-Δ(u(t))(x)

# u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
# u(t::Real) = x -> u(x,t)
# v(x) = t -> u(x,t)
# ∂tu(t) = x -> ForwardDiff.derivative(v(x),t)
# ∂tu(x,t) = ∂tu(t)(x)
# ∂t(::typeof(u)) = ∂tu
# f(t) = x -> ∂t(u)(t)(x)-Δ(u(t))(x)

p(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
p(t::Real) = x -> p(x,t)
q(x) = t -> p(x,t)
∂tp(t) = x -> ForwardDiff.derivative(q(x),t)
∂tp(x,t) = ∂tp(t)(x)
∂t(::typeof(p)) = ∂tp
g(t) = x -> ∂t(p)(t)(x)-Δ(p(t))(x)

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

order = 2

V0 = FESpace(
  reffe=:Lagrangian, order=order, valuetype=VectorValue{2,Float64},
  conformity=:H1, model=model, dirichlet_tags="boundary")
# V0 = FESpace(
#   reffe=:Lagrangian, order=order, valuetype=Float64,
#   conformity=:H1, model=model, dirichlet_tags="boundary")
U = TransientTrialFESpace(V0,u)

Q = FESpace(
  reffe=:Lagrangian, order=order, valuetype=Float64,
  conformity=:H1, model=model, dirichlet_tags="boundary")
P = TransientTrialFESpace(Q,p)

trian = Triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

#
a(u,v) = inner(∇(v),∇(u))
b(v,t) = inner(v,f(t))

X = MultiFieldFESpace([U,P])
Y = MultiFieldFESpace([V0,Q])

function res(t,x,xt,y)
  u,p = x
  ut,pt = xt
  v,q = y
  a(u,v) + inner(ut,v) - inner(v,f(t)) + a(p,q) + pt*q - q*g(t)
end

function jac(t,x,xt,dx,y)
  du1,du2 = dx
  v1,v2 = y
  a(du1,v1)+a(du2,v2)
end

function jac_t(t,x,xt,dxt,y)
  du1t,du2t = dxt
  v1,v2 = y
  inner(du1t,v1)+inner(du2t,v2)
end

function b(y)
  v,q = y
  0.0
  v*f(0.0) + q*g(0.0)
end

function mat(dx,y)
  du1,du2 = dx
  v1,v2 = y
  a(du1,v1)+a(du2,v2)
end

U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)
uh0 = interpolate_everywhere(U0,u(0.0))
ph0 = interpolate_everywhere(P0,p(0.0))
xh0 = Gridap.MultiField.MultiFieldFEFunction(X0,[uh0,ph0])

# t_Ω = AffineFETerm(mat,b,trian,quad)
# op = AffineFEOperator(X(0.0),Y,t_Ω)
# ls = LUSolver()
# solver = LinearFESolver(ls)
# uh = solve(solver,op)
# yh0 = Gridap.MultiField.MultiFieldFEFunction(X0,[uh0,ph0])
# algop = get_algebraic_operator(op)
# uh = solve!(yh0.free_values,ls,algop,nothing)

t_Ω = FETerm(res,jac,jac_t,trian,quad)
op = TransientFEOperator(X,Y,t_Ω)

t0 = 0.0
tF = 1.0
dt = 0.1

ls = LUSolver()
# using Gridap.Algebra: NewtonRaphsonSolver
# nls = NLSolver(ls;show_trace=true,method=:newton) #linesearch=BackTracking())
odes = ThetaMethod(ls,dt,θ)
solver = TransientFESolver(odes)

sol_t = solve(solver,op,xh0,t0,tF)

l2(w) = w*w


tol = 1.0e-6
_t_n = t0

result = Base.iterate(sol_t)

for (xh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  uh_tn = xh_tn.blocks[1]
  ph_tn = xh_tn.blocks[2]
  e = u(tn) - uh_tn
  el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
  e = p(tn) - ph_tn
  el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
  @test el2 < tol
end

# end #module
