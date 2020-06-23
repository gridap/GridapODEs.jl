# module DiffEqsWrapperTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator
using LineSearches: BackTracking

import Gridap: ∇
import GridapODEs.TransientFETools: ∂t

θ = 1.0

u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(t::Real) = x -> u(x,t)
∂tu = ∂t(u)

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
a(u,v) = ∇(v)⋅∇(u)
b(v,t) = v*f(t)
m(ut,v) = ut*v

#res(t,u,ut,v) = a(u,v) + ut*v - b(v,t)
#jac(t,u,ut,du,v) = a(du,v)
#jac_t(t,u,ut,dut,v) = dut*v

#t_Ω = FETerm(res,jac,jac_t,trian,quad)
t_Ω = TransientAffineFETerm(m,a,b,trian,quad)
op = TransientFEOperator(U,V0,t_Ω)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(U0,u(0.0))

ls = LUSolver()
using Gridap.Algebra: NewtonRaphsonSolver
nls = NLSolver(
        show_trace = true,
        method = :newton,
        linesearch = BackTracking(),
    )
# nls = NLSolver(ls;show_trace=true,method=:newton)
odes = ThetaMethod(nls,dt,θ)
solver = TransientFESolver(odes)
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
