module VectorHeatEquationTests

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

u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(t::Real) = x -> u(x,t)
v(x) = t -> u(x,t)
f(t) = x -> ∂t(u)(t)(x)-Δ(u(t))(x)

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

order = 2

V0 = FESpace(
  reffe=:Lagrangian, order=order, valuetype=Float64,
  conformity=:H1, model=model, dirichlet_tags="boundary")
U = TransientTrialFESpace(V0,u)

trian = Triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

#
a(u,v) = ∇(v)⊙∇(u)
b(v,t) = v⋅f(t)

X = TransientMultiFieldFESpace([U,U])
Y = MultiFieldFESpace([V0,V0])

_res(t,u,ut,v) = a(u,v) + ut⋅v - b(v,t)
_jac(t,u,ut,du,v) = a(du,v)
_jac_t(t,u,ut,dut,v) = dut⋅v

function res(t,x,xt,y)
  u1,u2 = x
  u1t,u2t = xt
  v1,v2 = y
  _res(t,u1,u1t,v1) + _res(t,u2,u2t,v2)
end

function jac(t,x,xt,dx,y)
  du1,du2 = dx
  v1,v2 = y
  a(du1,v1)+a(du2,v2)
end

function jac_t(t,x,xt,dxt,y)
  du1t,du2t = dxt
  v1,v2 = y
  du1t⋅v1+du2t⋅v2
end

t_Ω = FETerm(res,jac,jac_t,trian,quad)
op = TransientFEOperator(X,Y,t_Ω)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
X0 = X(0.0)
uh0 = interpolate_everywhere(U0,u(0.0))
xh0 = Gridap.MultiField.MultiFieldFEFunction(X0,[uh0,uh0])

ls = LUSolver()
# using Gridap.Algebra: NewtonRaphsonSolver
# nls = NLSolver(ls;show_trace=true,method=:newton) #linesearch=BackTracking())
odes = ThetaMethod(ls,dt,θ)
solver = TransientFESolver(odes)

sol_t = solve(solver,op,xh0,t0,tF)

l2(w) = w⋅w


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

end #module
