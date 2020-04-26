module HeatVectorEquationTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator

import Gridap: ∇
import GridapODEs.TransientFETools: ∂t

θ = 0.5

u(x,t) = VectorValue(x[1],x[2])*t
u(t::Real) = x -> u(x,t)
∂tu(t) = x -> VectorValue(x[1],x[2])
∂tu(x,t) = ∂tu(t)(x)
∂t(::typeof(u)) = ∂tu

f(t) = x -> ∂t(u)(t)(x)-Δ(u(t))(x)

domain = (0,1,0,1)
partition = (1,1)
model = CartesianDiscreteModel(domain,partition)

order = 2

V0 = FESpace(
  reffe=:Lagrangian, order=order, valuetype=VectorValue{2,Float64},
  conformity=:H1, model=model, dirichlet_tags="boundary")

U = TransientTrialFESpace(V0,u)

trian = Triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

a(u,v) = inner(∇(v),∇(u))
m(u,v) = inner(u,v)
b(v,t) = inner(v,f(t))

res(t,u,ut,v) = a(u,v) + m(ut,v) - b(v,t)
jac(t,u,ut,du,v) = a(du,v)
jac_t(t,u,ut,dut,v) = m(dut,v)

t_Ω = FETerm(res,jac,jac_t,trian,quad)
op = TransientFEOperator(U,V0,t_Ω)

t0 = 0.0
tF = 1.0
dt = 0.1


U0 = U(0.0)
uh0 = interpolate_everywhere(U0,u(0.0))

ls = LUSolver()
odes = ThetaMethod(ls,dt,θ)
solver = TransientFESolver(odes)

sol_t = solve(solver,op,uh0,t0,tF)

l2(w) = w*w

tol = 1.0e-6
_t_n = t0

rv, _ = Base.iterate(sol_t)
uh_tn, tn = rv
uh_tn.free_values

_t_n = t0
for (uh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  @test tn≈_t_n
  e = u(tn) - uh_tn
  el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
  @test el2 < tol
end

u0 = get_free_values(uh0)
uf = copy(get_free_values(uh0))

odeop = get_algebraic_operator(op)

ode_cache = allocate_cache(odeop)
vθ = similar(u0)
nl_cache = nothing

odes.θ == 0.0 ? dtθ = dt : dtθ = dt*odes.θ
tθ = t0+dtθ
ode_cache = update_cache!(ode_cache,odeop,tθ)

using GridapODEs.ODETools: ThetaMethodNonlinearOperator
nlop = ThetaMethodNonlinearOperator(odeop,tθ,dtθ,u0,ode_cache,vθ)

nl_cache = solve!(uf,odes.nls,nlop,nl_cache)
uf

K = nl_cache.A
h = nl_cache.b

# Steady version of the problem to extract the Laplacian and mass matrices
tf = tθ
Utf = U(tf)
fst(x) = f(tf)(x)
a(u,v) = inner(∇(v),∇(u))

function extract_matrix_vector(a,fst)
  btf(v) = inner(v,fst)
  t_Ω = AffineFETerm(a,btf,trian,quad)
  op = AffineFEOperator(Utf,V0,t_Ω)
  ls = LUSolver()
  solver = LinearFESolver(ls)
  uh = solve(solver,op)

  tol = 1.0e-6
  e = uh-u(tf)
  l2(e) = inner(e,e)
  l2e = sqrt(sum( integrate(l2(e),trian,quad) ))
  # @test l2e < tol

  Ast = op.op.matrix
  bst = op.op.vector

  @test uh.free_values ≈ Ast \ bst

  return Ast, bst
end

A,rhs = extract_matrix_vector(a,fst)

gst(x) = -1.0*u(0.0)(x)
m(u,v) = (1/(θ*dt))inner(u,v)

M,rhs2 = extract_matrix_vector(m,gst)

@test rhs + rhs2 ≈ h
@test A+M ≈ K
@test K \ h ≈ uf
@test K \ h ≈ uf

end #module
