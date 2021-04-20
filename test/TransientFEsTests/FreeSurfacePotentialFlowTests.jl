module FreeSurfacePotentialFlowTests

using Gridap
using Gridap.Geometry
using GridapODEs.TransientFETools
using GridapODEs.ODETools
using Test

# Parameters
L = 2*π
H = 1.0
n = 8
order = 2
g = 9.81
ξ = 0.1
λ = L/2
k = 2*π/L
h = L/n
ω = √(g*k*tanh(k*H))
t₀ = 0.0
tf = 2*π
Δt = h/(2*λ*ω)
θ = 0.5

# Exact solution
ϕₑ(x,t) = ω/k * ξ * (cosh(k*(x[2]))) / sinh(k*H) * sin(k*x[1] - ω*t)
ηₑ(x,t) = ξ * cos(k*x[1] - ω*t)
ϕₑ(t::Real) = x -> ϕₑ(x,t)
ηₑ(t::Real) = x -> ηₑ(x,t)

# Domain
domain = (0,L,0,H)
partition = (n,n)
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false))

# Boundaries
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"free_surface",[3,4,6])
bgface_to_mask = get_face_mask(labels,[3,4,6],1)
model_Γ =BoundaryDiscreteModel(Polytope{1},model,bgface_to_mask)

# Triangulation
Ω = Triangulation(model)
Γ = Triangulation(model_Γ)
dΩ = Measure(Ω,2*order)
dΓ = Measure(Γ,2*order)

# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,conformity=:H1)
V_Γ = TestFESpace(model_Γ,reffe,conformity=:H1)
U = TransientTrialFESpace(V)
U_Γ = TransientTrialFESpace(V_Γ)
X = TransientMultiFieldFESpace([U,U_Γ])
Y = MultiFieldFESpace([V,V_Γ])

# Weak form
α = 2/Δt

# Optimal transient FE Operator:
m((ϕt,ηt),(w,v)) = ∫( 0.5*(α/g*(w*ϕt) + v*ϕt) - (w*ηt) )dΓ
a((ϕ,η),(w,v)) = ∫( ∇(ϕ)⋅∇(w) )dΩ + ∫( 0.5*(α*(w*η) + g*v*η) )dΓ
b((w,v)) = ∫( 0.0*w )dΓ
# op = TransientConstantFEOperator(m,a,b,X,Y)

# TransientFEOperator exploiting automatic differentiation (testing purposes)
res(t,x,y) = m(∂t(x),y) + a(x,y) - b(y)
jac(t,x,dx,y) = a(dx,y)
jac_t(t,x,dxt,y) = m(dxt,y)
op = TransientFEOperator(res,jac,jac_t,X,Y)

# Solver
ls = LUSolver()
odes = ThetaMethod(ls,Δt,θ)
solver = TransientFESolver(odes)

# Initial solution
x₀ = interpolate_everywhere([ϕₑ(0.0),ηₑ(0.0)],X(0.0))

# Solution
sol_t = solve(solver,op,x₀,t₀,tf)

# Post-process
l2_Ω(w) = √(∑(∫(w*w)dΩ))
l2_Γ(v) = √(∑(∫(v*v)dΓ))
E_kin(w) = 0.5*∑( ∫(∇(w)⋅∇(w))dΩ )
E_pot(v) = g*0.5*∑( ∫(v*v)dΓ )
Eₑ = 0.5*g*ξ^2*L

tol = 1.0e-2
for ((ϕn,ηn),tn) in sol_t
  E = E_kin(ϕn) + E_pot(ηn)
  error_ϕ = l2_Ω(ϕn-ϕₑ(tn))
  error_η = l2_Γ(ηn-ηₑ(tn))
  @test abs(E/Eₑ-1.0) <= tol
  @test error_ϕ <= tol
  @test error_η <= tol
end

end
