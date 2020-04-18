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

θ = 0.5

# Analytical functions
# u(x,t) = (x[1]+x[2])*t
# u(x,t) = (2*x[1]+x[2])*t
u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t

# @santiagobadia : @fverdugo, take a look at this, we could probably start using
# it in the tutorials, etc, creating functions that return gradients, laplacians,
# etc
u(t::Real) = x -> u(x,t)

∇u(x,t) = VectorValue(ForwardDiff.gradient(u(t),x)...)
∇u(x::Gridap.TensorValues.MultiValue,t) = ∇u(x.array,t)
∇u(t::Real) = x -> ∇u(x,t)

_∇u(x,t) = ForwardDiff.gradient(u(t),x)
_∇u(t::Real) = x -> ForwardDiff.gradient(u(t),x)
Δu(x,t) = tr(ForwardDiff.jacobian(_∇u(t),x))
Δu(x::Gridap.TensorValues.MultiValue,t) = Δu(x.array,t)
Δu(t::Real) = x -> Δu(x,t)

v(x) = t -> u(x,t)
∂tu(t) = x -> ForwardDiff.derivative(v(x),t)
∂tu(x,t) = ∂tu(t)(x)

∂t(::typeof(u)) = ∂tu
∇(::typeof(u)) = ∇u
#

f(t) = x -> ∂tu(x,t)-Δu(x,t)

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

order = 2

V0 = FESpace(
  reffe=:Lagrangian, order=order, valuetype=Float64,
  conformity=:H1, model=model, dirichlet_tags="boundary")

U = TransientTrialFESpace(V0,u)

X = MultiFieldFESpace([U,U])
Y = MultiFieldFESpace([V0,V0])

U.Ud0

MultiFieldFESpace([fesp.Ud0 for fesp in X.spaces])

trian = Triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

#
a(u,v) = ∇(v)*∇(u)
b(v,t) = v*f(t)

_res(t,u,ut,v) = a(u,v) + ut*v - b(v,t)
_jac(t,u,ut,du,v) = a(du,v)
_jac_t(t,u,ut,dut,v) = dut*v


function res(t,x,xt,y)
  u1,u2 = x
  u1t,u2t = xt
  v1,v2 = y
  a(u1,v1) + u1t*v1 - b(v1,t) + a(u2,v2) + u2t*v2 - b(v2,t)
end

function jac(t,x,xt,dx,y)
  du1,du2 = dx
  v1,v2 = y
  a(du1,v1)+a(du2,v2)
end

function jac_t(t,x,xt,dxt,y)
  du1t,du2t = dxt
  v1,v2 = y
  du1t*v1+du2t*v2
end

t_Ω = FETerm(res,jac,jac_t,trian,quad)
op = TransientFEOperator(X,Y,t_Ω)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(U0,u(0.0))

X0 = X(0.0)
uh0 = interpolate_everywhere(X0,u(0.0))

MultiFi

ls = LUSolver()
using Gridap.Algebra: NewtonRaphsonSolver
nls = NLSolver(ls;show_trace=true,method=:newton) #linesearch=BackTracking())
odes = ThetaMethod(ls,dt,θ)
solver = TransientFESolver(odes) # Return a specialization of TransientFESolver

sol_t = solve(solver,op,uh0,t0,tF)

l2(w) = w*w


tol = 1.0e-6
_t_n = t0

for (uh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  e = u(tn) - uh_tn
  el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
  @test el2 < tol
end

# end #module
