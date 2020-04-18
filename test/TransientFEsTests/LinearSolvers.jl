
using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator

import Gridap: ∇
import GridapODEs.TransientFETools: ∂t

θ = 1.0

# Analytical functions
u(x,t) = (x[1]+x[2])*t
# u(x,t) = (2*x[1]+x[2])*t
# u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t

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

trian = Triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

a(u,v) = ∇(v)*∇(u)
b(v,t) = v*f(t)

res(t,u,ut,v) = a(u,v) + ut*v - b(v,t)
jac(t,u,ut,du,v) = a(du,v)
jac_t(t,u,ut,dut,v) = dut*v

t_Ω = FETerm(res,jac,jac_t,trian,quad)
op = TransientFEOperator(U,V0,t_Ω)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(U0,u(0.0))

ls = LUSolver()

odes = ThetaMethod(ls,dt,θ)

solver = TransientFESolver(odes) # Return a specialization of TransientFESolver

# solver = odes
# u0 = get_free_values(uh0)
# odeop = get_algebraic_operator(op)
# ode_cache = allocate_cache(odeop)
# uf = copy(u0)
# tf = t0+dt
# cache = nothing
#
# ode_cache = allocate_cache(odeop)
# vθ = similar(u0)
# ls_cache = nothing
#
# dt = solver.dt
# tf = t0+dt
# solver.θ == 0.0 ? dtθ = dt : dtθ = dt*solver.θ
# tθ = t0+dtθ
# ode_cache = update_cache!(ode_cache,odeop,tf)
#
# X = GridapODEs.ODETools
# lop = X.ThetaMethodNonlinearOperator(odeop,tθ,dtθ,u0,ode_cache,vθ)
# A = X.allocate_jacobian(lop,lop.u0)
# z = similar(lop.u0)
# z .= 0.0
# bb = X.allocate_residual(lop,z)
# # A, b = ls_cache
# jacobian!(A,lop,lop.u0)
# X.residual!(bb,lop,z)
# X.AffineOperator(A,-bb)
#
# lop = GridapODEs.ODETools.ThetaMethodOperator(solver,odeop,tθ,dtθ,u0,ode_cache,vθ,ls_cache)

#
# u0 = -J1\b1
# t0 = dt
# tf = t0+dt
#
# nlop = BackwardEulerNonlinearOperator(odeop,tf,dt,u0,ode_cache)
# nlop = BackwardEulerNonlinearOperator(odeop,tf,dt,u0,ode_cache)
#
# x = copy(nlop.u0)
# b1 = allocate_residual(nlop,x)
# residual!(b1,nlop,x)
# J1 = allocate_jacobian(nlop,x)
# jacobian!(J1,nlop,x)


# t0 = 0.0
# tf = t0+dt

sol_t = solve(solver,op,uh0,t0,tF)


l2(w) = w*w


tol = 1.0e-6
_t_n = t0

for (uh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  # @test tn≈_t_n
  e = u(tn) - uh_tn
  @show tn
  @show get_free_values(uh_tn)
  @show get_dirichlet_values(uh_tn)
  el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
  @test el2 < tol
end

# end #module
