using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using GridapTimeStepper.ODETools
using GridapTimeStepper.TransientFETools
using Gridap.FESpaces: get_algebraic_operator

import Gridap: ∇
import GridapTimeStepper.TransientFETools: ∂t

θ = 0.4

# Analytical functions
u(x,t) = (x[1]+x[2])*t
u(x,t) = (2*x[1]+x[2])*t
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

trian = Triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

#
a(u,v) = ∇(v)*∇(u)
b(v,t) = v*f(t)


# g(x) = ∂tu(x,2.0)
# g(x) = -Δu(x,1.0)
# b(v,t) = v*g

# Steady version of the problem to extract the Laplacian and mass matrices

# U0 = U(1.0)
# _res_lap(u,v) = a(u,v) - b(v,0.0)
# _jac_lap(u,du,v) = a(du,v)
# _mas(u,du,v) = du*v
# # _t_lap = FETerm(_res_lap,_jac_lap,trian,quad)
# _t_lap = FETerm(_res_lap,_mas,trian,quad)
# _op = FEOperator(U0,V0,_t_lap)
# uh = interpolate_everywhere(U0,0.0)#1.0)
# using Gridap.FESpaces: allocate_residual, allocate_jacobian
# _b = allocate_residual(_op,uh)
# # Juno.@enter allocate_residual(_op,uh)
# _L = allocate_jacobian(_op,uh)
# using Gridap.FESpaces: residual!, jacobian!
# residual!(_b,_op,uh)
# jacobian!(_L,_op,uh)
# @show _b
# @show _L
# @show _L \ _b
#

# b0(v) = 0.0*b(v,0.0)
# t_Ω = AffineFETerm(a,b0,trian,quad)
# _op = AffineFEOperator(U0,V0,t_Ω)
# ls = LUSolver()
# solver = LinearFESolver(ls)
# uh = solve(solver,_op)
# writevtk(trian,"results",cellfields=["uh"=>uh])
# @show get_free_values(uh)
# @show uh
# end steady


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
# u0 = get_free_values(uh0)
#
# u0 = get_free_values(uh0)
# odeop = get_algebraic_operator(op)
# ode_cache = allocate_cache(odeop)
# uf = copy(u0)
# tf = t0+dt
# update_cache!(ode_cache,odeop,tf)
# using GridapTimeStepper.ODETools: BackwardEulerNonlinearOperator
# nlop = BackwardEulerNonlinearOperator(odeop,tf,dt,u0,ode_cache)
# using GridapTimeStepper.ODETools: BackwardEulerNonlinearOperator
# nlop = BackwardEulerNonlinearOperator(odeop,tf,dt,u0,ode_cache)
#
# x = copy(nlop.u0)
# b1 = allocate_residual(nlop,x)
# residual!(b1,nlop,x)
# J1 = allocate_jacobian(nlop,x)
# jacobian!(J1,nlop,x)
#
ls = LUSolver()
using Gridap.Algebra: NewtonRaphsonSolver
nls = NLSolver(ls;show_trace=true,method=:newton) #linesearch=BackTracking())
odes = ThetaMethod(nls,dt,θ)
# odes = BackwardEuler(nls,dt)
solver = TransientFESolver(odes) # Return a specialization of TransientFESolver
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

# _t_n = t0
# result, state = Base.iterate(sol_t)
# uh_tn = result[1]
# @show get_free_values(uh)
# @show get_dirichlet_values(uh)
# tn = result[2]
# u(Point(0.5,0.5),0.1)
# e = u(tn) - uh_tn
# el2 = sqrt(sum( integrate(l2(e),trian,quad) ))

# un, tn = Base.iterate(sol_t)

# result, state = Base.iterate(sol_t)
# get_free_values(result[1])
# result, state = Base.iterate(sol_t,state)
# get_free_values(result[1])
# result[2]

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
