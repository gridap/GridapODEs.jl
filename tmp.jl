using Gridap
using Gridap.Arrays
using Gridap.Geometry
using GridapODEs.TransientFETools
using GridapODEs.ODETools
using Test
using LineSearches: BackTracking
using ForwardDiff
import GridapODEs.TransientFETools: ∂t

# Problem setting
println("Setting analytical fluid problem parameters")
# Solid properties
E_s = 1.0
ν_s = 0.4
ρ_s = 1.0
# Fluid properties
ρ_f = 1.0
μ_f = 1.0
γ_f = 1.0
# Mesh properties
n_m = 3
E_m = 1.0
ν_m = 0.1
α_m = 1.0e-5
# Time stepping
t0 = 0.0
tf = 0.5
dt = 0.1
θ  = 0.5

# Define BC functions
println("Defining Boundary conditions")
#u(x, t) = 1.0e-2 * VectorValue( x[1]^2*x[2] , -x[1]*x[2]^2) * t
#v(x, t) = 1.0e-2 * VectorValue( x[1]^2*x[2] , -x[1]*x[2]^2)
u(x, t) = 1.0e-2 * VectorValue( 1.0 , -1.0) * t
v(x, t) = 1.0e-2 * VectorValue( 1.0 , -1.0)
u(t::Real) = x -> u(x,t)
v(t::Real) = x -> v(x,t)
∂tu(t) = x -> VectorValue(ForwardDiff.derivative(t -> get_array(u(x,t)),t))
∂tv(t) = x -> VectorValue(ForwardDiff.derivative(t -> get_array(v(x,t)),t))
∂tu(x,t) = ∂tu(t)(x)
∂tv(x,t) = ∂tv(t)(x)
T_u = typeof(u)
T_v = typeof(v)
@eval ∂t(::$T_u) = $∂tu
@eval ∂t(::$T_v) = $∂tv
p(x,t) = 0.0
p(t::Real) = x -> p(x,t)
bconds = (
  # Tags
  FSI_Vu_tags = ["boundary"],
  FSI_Vv_tags = ["boundary"],
  ST_Vu_tags = ["boundary","interface"],
  ST_Vv_tags = ["boundary","interface"],
  # Values,
  FSI_Vu_values = [u],
  FSI_Vv_values = [v],
  ST_Vu_values = [u(0.0),u(0.0)],
  ST_Vv_values = [v(0.0),v(0.0)],
)

function lame_parameters(E,ν)
  λ = (E*ν)/((1+ν)*(1-2*ν))
  μ = E/(2*(1+ν))
  (λ, μ)
end

# Define Forcing terms
I = TensorValue( 1.0, 0.0, 0.0, 1.0 )
F(t) = x -> ∇(u(t))(x) + I
J(t) = x -> det(F(t))(x)
E(t) = x -> 0.5 * ((F(t)(x)')⋅F(t)(x) - I)
(λ_s, μ_s) = lame_parameters(E_s,ν_s)
(λ_m, μ_m) = lame_parameters(E_m,ν_m)
S_SV(t) = x -> 2*μ_s*E(t)(x) + λ_s*tr(E(t)(x))*I
fv_ST_Ωf(t) = x -> - μ_f*Δ(v(t))(x) + ∇(p(t))(x)
fu_Ωf(t) = x -> - α_m * Δ(u(t))(x)
fv_Ωf(t) = x -> ρ_f * ∂t(v)(t)(x) - μ_f * Δ(v(t))(x) + ∇(p(t))(x) + ρ_f*( (∇(v(t))(x)')⋅(v(t)(x) - ∂t(u)(t)(x)) )
fp_Ωf(t) = x -> (∇⋅v(t))(x)
fu_Ωs(t) = x -> ∂t(u)(t)(x) - v(t)(x)
fv_Ωs(t) = x -> ρ_s * ∂t(v)(t)(x) #- (∇⋅(F(t)⋅S_SV(t)))(x)  # Divergence of a a doted function not working yet...

# Discrete model
println("Defining discrete model")
domain = (-1,1,-1,1)
partition = (n_m,n_m)
model = CartesianDiscreteModel(domain,partition)
trian = Triangulation(model)
R = 0.5
xc = 0.0
yc = 0.0
function is_in(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  d = (x[1]-xc)^2 + (x[2]-yc)^2 - R^2
  d < 1.0e-8
end
oldcell_to_coods = get_cell_coordinates(trian)
oldcell_to_is_in = collect1d(lazy_map(is_in,oldcell_to_coods))
incell_to_cell = findall(oldcell_to_is_in)
outcell_to_cell = findall(collect(Bool, .! oldcell_to_is_in))
model_solid = DiscreteModel(model,incell_to_cell)
model_fluid = DiscreteModel(model,outcell_to_cell)

# Build fluid-solid interface labelling
println("Defining Fluid-Solid interface")
labeling = get_face_labeling(model_fluid)
new_entity = num_entities(labeling) + 1
topo = get_grid_topology(model_fluid)
D = num_cell_dims(model_fluid)
for d in 0:D-1
  fluid_boundary_mask = collect(Bool,get_isboundary_face(topo,d))
  fluid_outer_boundary_mask = get_face_mask(labeling,"boundary",d)
  fluid_interface_mask = collect(Bool,fluid_boundary_mask .!= fluid_outer_boundary_mask)
  dface_list = findall(fluid_interface_mask)
  for face in dface_list
    labeling.d_to_dface_to_entity[d+1][face] = new_entity
  end
end
add_tag!(labeling,"interface",[new_entity])

# Triangulations
println("Defining triangulations")
Ω = Triangulation(model)
Ω_s = Triangulation(model_solid)
Ω_f = Triangulation(model_fluid)
Γ_i = BoundaryTriangulation(model_fluid,tags="interface")

# Quadratures
println("Defining quadratures")
order = 2
degree = 2*order
bdegree = 2*order
dΩ  = Measure(Ω,degree)
dΩs = Measure(Ω_s,degree)
dΩf = Measure(Ω_f,degree)
dΓi = Measure(Γ_i,bdegree)

# Test FE Spaces
println("Defining FE spaces")
# Reference FE
reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
reffeₚ = ReferenceFE(lagrangian,Float64,order-1)
# Test FE Spaces
Vu_FSI = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=bconds[:FSI_Vu_tags])
Vv_FSI = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=bconds[:FSI_Vv_tags])
Vu_ST = TestFESpace(model_fluid, reffeᵤ, conformity=:H1, dirichlet_tags=bconds[:ST_Vu_tags])
Vv_ST = TestFESpace(model_fluid, reffeᵤ, conformity=:H1, dirichlet_tags=bconds[:ST_Vv_tags])
Q = TestFESpace(model_fluid, reffeₚ, conformity=:C0, constraint=:zeromean)

# Trial FE Spaces
Uu_ST = TrialFESpace(Vu_ST,bconds[:ST_Vu_values])
Uv_ST = TrialFESpace(Vv_ST,bconds[:ST_Vv_values])
Uu_FSI = TransientTrialFESpace(Vu_FSI,bconds[:FSI_Vu_values])
Uv_FSI = TransientTrialFESpace(Vv_FSI,bconds[:FSI_Vv_values])
P = TrialFESpace(Q)

# Multifield FE Spaces
Y_ST = MultiFieldFESpace([Vu_ST,Vv_ST,Q])
X_ST = MultiFieldFESpace([Uu_ST,Uv_ST,P])
Y_FSI = MultiFieldFESpace([Vu_FSI,Vv_FSI,Q])
X_FSI = TransientMultiFieldFESpace([Uu_FSI,Uv_FSI,P])

# Stokes problem for initial solution
println("Defining Stokes operator")
a((u,v,p),(ϕ,φ,q)) = ∫( ϕ⋅u  + ε(φ) ⊙ (2.0*μ_f*ε(v)) - (∇⋅φ) * p + q * (∇⋅v) )dΩf
l((ϕ,φ,q)) = ∫( φ⋅fv_ST_Ωf(0.0) )dΩf
res(x,y) = a(x,y) - l(y)
jac(x,dx,y) = a(dx,y)
op_ST = FEOperator(res,jac,X_ST,Y_ST)

# Map operations
F(∇u) = ∇u + I
J(∇u) = det(F(∇u))
Finv(∇u) = inv(F(∇u))
FinvT(∇u) = (Finv(∇u)')
E(∇u) = 0.5 * ((F(∇u)')⋅F(∇u) - I)

# Map operations derivatives
dF(∇du) = ∇du
dJ(∇u,∇du) = J(∇u)*tr(Finv(∇u)⋅dF(∇du))
dFinv(∇u,∇du) = -Finv(∇u) ⋅ dF(∇du) ⋅ Finv(∇u)
dFinvT(∇u,∇du) = (dFinv(∇u,∇du)')
dE(∇u,∇du) = 0.5 * ((dF(∇du)')⋅F(∇u) + (F(∇u)')⋅dF(∇du))

# Fluid constitutive laws
# Cauchy stress
σᵥ_Ωf(μ,u) = 2.0*μ*ε(u)
σᵥ_Ωf(μ) = (∇v,Finv) -> μ*(∇v⋅Finv + (Finv')⋅(∇v'))
# First Piola-Kirchhoff stress
Pᵥ_Ωf(μ,u,v) = (J∘∇(u)) * (σᵥ_Ωf(μ)∘(∇(v),Finv∘∇(u)))' ⋅ (FinvT∘∇(u))
function Pₚ_Ωf(u,p)
  typeof(- (J∘∇(u)) * p * tr((FinvT∘∇(u))))
  - (J∘∇(u)) * p * tr((FinvT∘∇(u)))
end
# First Piola-Kirchhoff stress Jacobian
dPᵥ_Ωf_du(μ,u,du,v) = (dJ∘(∇(u),∇(du))) * (σᵥ_Ωf(μ)∘(∇(v),(Finv∘∇(u))))' ⋅ (FinvT∘∇(u)) +
                      (J∘∇(u)) * (σᵥ_Ωf(μ)∘(∇(v),(dFinv∘(∇(u),∇(du)))))' ⋅ (FinvT∘∇(u)) +
                      (J∘∇(u)) * (σᵥ_Ωf(μ)∘(∇(v),(Finv∘∇(u))))' ⋅ (dFinvT∘(∇(u),∇(du)))
dPₚ_Ωf_du(u,du,p) = - p * ( (dJ∘(∇(u),∇(du))) * tr((FinvT∘∇(u))) +
                           (J∘∇(u)) * tr((dFinvT∘(∇(u),∇(du)))) )
dPᵥ_Ωf_dv(μ,u,dv) = Pᵥ_Ωf(μ,u,dv)
dPₚ_Ωf_dp(u,dp) = Pₚ_Ωf(u,dp)
# Convective term
conv(c,∇v) = (∇v') ⋅ c
dconv(dc,∇dv,c,∇v) = conv(c,∇dv) + conv(dc,∇v)

# Solid constitutive laws
# Second Piola-Kirchhoff stress (Saint-Venant)
Sₛᵥ_Ωs(λ,μ,u) = 2*μ*(E∘∇(u)) + λ*tr((E∘∇(u)))*(I∘∇(u))
dSₛᵥ_Ωs_du(λ,μ,u,du) = 2*μ*(dE∘(∇(u),∇(du))) + λ*tr((dE∘(∇(u),∇(du))))*(I∘∇(u))
# First Piola-Kirchhoff stress (Saint-Venant)
Pₛᵥ_Ωs(λ,μ,u) = (F∘∇(u)) ⋅ Sₛᵥ_Ωs(λ,μ,u)
dPₛᵥ_Ωs_du(λ,μ,u,du) = (dF∘∇(du)) ⋅ Sₛᵥ_Ωs(λ,μ,u) + (F∘∇(u)) ⋅ dSₛᵥ_Ωs_du(λ,μ,u,du)

# # FSI problem
println("Defining FSI operator")
n_Γi = get_normal_vector(Γ_i)
a_f((u,v,p),(ϕ,φ,q)) = ∫( φ ⋅ ( (J∘∇(u)) * ρ_f * ∂t(v) ) )dΩ 
da_f((u,v,p),(du,dv,dp),(ϕ,φ,q)) =  ∫( φ ⋅ ( (dJ∘(∇(u),∇(du))) * ρ_f * ∂t(v) ) )dΩ 
da_t_f((u,v,p),(dut,dvt,dpt),(ϕ,φ,q)) = ∫( φ ⋅ ( (J∘∇(u)) * ρ_f * dvt ) )dΩ
# da_t_f((u,v,p),(dut,dvt,dpt),(ϕ,φ,q)) = ∫( φ ⋅ ( ρ_f * dvt) )dΩf
a_fsi((u,v,p),(ϕ,φ,q)) = ∫( ϕ⋅u + φ⋅v + q*p )dΩ
l_f(t,(ϕ,φ,q)) = ∫( ϕ⋅fu_Ωf(t) )dΩ
l_fsi(t,(ϕ,φ,q)) = ∫( φ⋅fv_Ωf(t) )dΩ
res(t,x,y) = a_f(x,y) + a_fsi(x,y) + a_fsi(∂t(x),y) - l_f(t,y) - l_fsi(t,y)
jac(t,x,dx,y) = da_f(x,dx,y) + a_fsi(dx,y)
jac_t(t,x,dxt,y) = da_t_f(x,dxt,y) + a_fsi(dxt,y)
op_FSI = TransientFEOperator(res,jac,jac_t,X_FSI,Y_FSI)
#op_FSI = TransientFEOperator(res,X_FSI,Y_FSI)
# a_f((u,v,p),(ut,vt,pt),(ϕ,φ,q)) =   
#   ∫( α_m * (∇(ϕ) ⊙ ∇(u)) )dΩf +
#   ∫( φ ⋅ ( (J∘∇(u)) * ρ_f * vt ) )dΩf 
# da_f((u,v,p),(ut,vt,pt),(du,dv,dp),(ϕ,φ,q)) =   
#   ∫( α_m * (∇(ϕ) ⊙ ∇(du)) )dΩf +
#   ∫( φ ⋅ ( (dJ∘(∇(u),∇(du))) * ρ_f * vt ) )dΩf 
# da_t_f((u,v,p),(dut,dvt,dpt),(ϕ,φ,q)) = ∫( φ ⋅ ( (J∘∇(u)) * ρ_f * dvt ) )dΩf
# a_fsi((u,v,p),(ϕ,φ,q)) = ∫( ϕ⋅u + φ⋅v + q*p )dΩf + ∫( ϕ⋅u + φ⋅v + q*p )dΩs
# l_f(t,(ϕ,φ,q)) = ∫( ϕ⋅fu_Ωf(t) )dΩf
# l_fsi(t,(ϕ,φ,q)) = ∫( φ⋅fv_Ωf(t) )dΩf
# res(t,(x,xt),y) = a_f(x,xt,y) + a_fsi(x,y) + a_fsi(xt,y) - l_f(t,y) - l_fsi(t,y)
# jac(t,(x,xt),dx,y) = da_f(x,xt,dx,y) + a_fsi(dx,y)
# jac_t(t,(x,xt),dxt,y) = da_t_f(x,dxt,y) + a_fsi(dxt,y)
# op_FSI = TransientFEOperator(res,jac,jac_t,X_FSI,Y_FSI)

# Setup output files
folderName = "fsi-results"
fileName = "fields"
if !isdir(folderName)
  mkdir(folderName)
end
filePath = join([folderName, fileName], "/")

# Solve Stokes problem
println("Defining Stokes solver")
xh = solve(op_ST)
writevtk(Ω_f, filePath, cellfields = ["uh"=>xh[1], "vh"=>xh[2], "ph"=>xh[3]])

# Compute Stokes solution L2-norm
l2(w) = w⋅w
eu_ST = u(0.0) - xh[1]
ev_ST = v(0.0) - xh[2]
ep_ST = p(0.0) - xh[3]
eul2_ST = sqrt(∑( ∫(l2(eu_ST))dΩf ))
evl2_ST = sqrt(∑( ∫(l2(ev_ST))dΩf ))
epl2_ST = sqrt(∑( ∫(l2(ep_ST))dΩf ))
println("Stokes L2-norm u: ", eul2_ST)
println("Stokes L2-norm v: ", evl2_ST)
println("Stokes L2-norm p: ", epl2_ST)
@test eul2_ST < 1.0e-10
@test evl2_ST < 1.0e-10
@test epl2_ST < 1.0e-10

# Solve FSI problem
println("Defining FSI solver")
xh0 = interpolate_everywhere([u(0.0),v(0.0),p(0.0)],X_FSI(0.0))
nls = NLSolver(
show_trace = true,
method = :newton,
linesearch = BackTracking(),
ftol = 1.0e-10,
iterations = 50
)
odes =  ThetaMethod(nls, dt, 0.5)
solver = TransientFESolver(odes)
xht = solve(solver, op_FSI, xh0, t0, tf)

for (i,(xh, t)) in enumerate(xht)
  println("STEP: $i, TIME: $t")
  println("============================")

  # Compute errors
  eu = u(t) - xh[1]
  ev = v(t) - xh[2]
  ep = p(t) - xh[3]
  eul2 = sqrt(∑( ∫(l2(eu))dΩ ))
  evl2 = sqrt(∑( ∫(l2(ev))dΩ ))
  epl2 = sqrt(∑( ∫(l2(ep))dΩ ))

  # Test checks
  println("FSI L2-norm u: ", eul2)
  println("FSI L2-norm v: ", evl2)
  println("FSI L2-norm p: ", epl2)
end