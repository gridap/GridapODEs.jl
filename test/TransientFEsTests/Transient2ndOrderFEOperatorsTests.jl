module Transient2nOrderFEOperatorsTests

using Gridap
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using Test

# Analytical functions
u(x,t) = (1.0-x[1])*x[1]*(t^2+3.0)
u(t::Real) = x -> u(x,t)
v(t::Real) = ∂t(u)(t)
a(t::Real) = ∂tt(u)(t)
f(t) = x -> ∂tt(u)(x,t) + ∂t(u)(x,t) - Δ(u(t))(x)

domain = (0,1)
partition = (2,)
model = CartesianDiscreteModel(domain,partition)

order = 2
reffe = ReferenceFE(lagrangian,Float64,order)
V0 = FESpace(
  model,
  reffe,
  conformity=:H1,
  dirichlet_tags="boundary")
U = TransientTrialFESpace(V0,u)

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)

m(utt,v) = ∫(v*utt)dΩ
c(ut,v) = ∫(v*ut)dΩ
a(u,v) = ∫(∇(v)⊙∇(u))dΩ
b(v,t) = ∫(v*f(t))dΩ

res(t,u,v) = m(∂tt(u),v) + c(∂t(u),v) + a(u,v) - b(v,t)
jac(t,u,du,v) = a(du,v)
jac_t(t,u,dut,v) = c(dut,v)
jac_tt(t,u,dutt,v) = m(dutt,v)

op = TransientFEOperator(res,jac,jac_t,jac_tt,U,V0)

t0 = 0.0
tF = 1.0
dt = 0.1
γ = 0.5
β = 0.25

U0 = U(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)
vh0 = interpolate_everywhere(v(0.0),U0)
ah0 = interpolate_everywhere(a(0.0),U0)

ls = LUSolver()
odes = Newmark(ls,dt,γ,β)
solver = TransientFESolver(odes)

sol_t = solve(solver,op,(uh0,vh0,ah0),t0,tF)

l2(w) = w*w

tol = 1.0e-6
_t_n = t0

for (uh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  @test tn≈_t_n
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
end

end
