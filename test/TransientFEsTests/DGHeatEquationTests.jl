module DGHeatEquationTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator

# import Gridap: ∇
# import GridapODEs.TransientFETools: ∂t

θ = 0.2

# Analytical functions
# u(x,t) = (x[1]+x[2])*t
# u(x,t) = (2*x[1]+x[2])*t
u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(t::Real) = x -> u(x,t)
# v(x) = t -> u(x,t)
# ∂tu(t) = x -> ForwardDiff.derivative(v(x),t)
# ∂tu(x,t) = ∂tu(t)(x)
# ∂t(::typeof(u)) = ∂tu
f(t) = x -> ∂t(u)(x,t)-Δ(u(t))(x)

L= 1.0
n = 2
domain = (0,L,0,L)
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

order = 2

V0 = FESpace(
  reffe=:Lagrangian, order=order, valuetype=Float64,
  conformity=:L2, model=model)
U = TransientTrialFESpace(V0)

# trian = Triangulation(model)
# degree = 2*order
# quad = CellQuadrature(trian,degree)
#
# #
# a(u,v) = ∇(v)⋅∇(u)
# b(v,t) = v*f(t)
#
# res(t,u,ut,v) = a(u,v) + ut*v - b(v,t)
# jac(t,u,ut,du,v) = a(du,v)
# jac_t(t,u,ut,dut,v) = dut*v
#
# t_Ω = FETerm(res,jac,jac_t,trian,quad)
#
# # neumanntags = [7,8]
# btrian = BoundaryTriangulation(model)
# # btrian = BoundaryTriangulation(model,neumanntags)
# bquad = CellQuadrature(btrian,degree)
# nb = get_normal_vector(btrian)
#
h = 1.0 / n
γ = order*(order+1)
# a_∂Ω(u,v) = (γ/h)*v*u - v*(∇(u)⋅nb) - (∇(v)⋅nb)*u
# b_∂Ω(v,t) = (γ/h)*v*u(t) - (∇(v)⋅nb)*u(t)
# # b_∂Ω(v,t) = v*(∇(u(t))⋅nb)
#
# res_∂Ω(t,u,ut,v) = a_∂Ω(u,v) - b_∂Ω(v,t)
# jac_∂Ω(t,u,ut,du,v) = a_∂Ω(du,v)
# jac_t_∂Ω(t,u,ut,dut,v) = dut*v*0.0
#
# # t_∂Ω = AffineFETerm(a_∂Ω,b_∂Ω,btrian,bquad)
# t_∂Ω = FETerm(res_∂Ω,jac_∂Ω,jac_t_∂Ω,btrian,bquad)


degree = 2*order
strian = SkeletonTriangulation(model)
squad = CellQuadrature(strian,degree)
ns = get_normal_vector(strian)

a_Γ(u,v) = (γ/h)*jump(v*ns)⊙jump(u*ns) - jump(v*ns)⊙mean(∇(u)) - mean(∇(v))⊙jump(u*ns)

res_Γ(t,u,ut,v) = a_Γ(u,v)
jac_Γ(t,u,ut,du,v) = a_Γ(du,v)
jac_t_Γ(t,u,ut,dut,v) = dut*v*0.0

t_Γ = FETerm(res_Γ,jac_Γ,jac_t_Γ,strian,squad)

# op = TransientFEOperator(U,V0,t_Ω,t_∂Ω,t_Γ)
op = TransientFEOperator(U,V0,t_Γ)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(U0,u(0.0))

ls = LUSolver()
using Gridap.Algebra: NewtonRaphsonSolver
# nls = NLSolver(ls;show_trace=true,method=:newton) #linesearch=BackTracking())
odes = ThetaMethod(ls,dt,θ)
solver = TransientFESolver(odes)

sol_t = solve(solver,op,uh0,t0,tF)

# Juno.@enter Base.iterate(sol_t)

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

end #module
