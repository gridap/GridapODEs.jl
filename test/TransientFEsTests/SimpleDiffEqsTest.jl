module GridapFESolver

using Test
using Gridap
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator
using GridapODEs

export fe_problem
export solve_ode_gridap
export diffeq_wrappers

# FE problem (heat eq) using Gridap
function fe_problem(u, n)

  f(t) = x -> ∂t(u)(x, t) - Δ(u(t))(x)

  domain = (0, 1, 0, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition)

  order = 1

  V0 = FESpace(
    reffe = :Lagrangian,
    order = order,
    valuetype = Float64,
    conformity = :H1,
    model = model,
    dirichlet_tags = "boundary",
  )
  U = TransientTrialFESpace(V0, u)

  trian = Triangulation(model)
  degree = 2 * order
  quad = CellQuadrature(trian, degree)

  a(u, v) = ∇(v) ⋅ ∇(u)
  b(v, t) = v * f(t)

  res(t, u, ut, v) = a(u, v) + ut * v - b(v, t)
  jac(t, u, ut, du, v) = a(du, v)
  jac_t(t, u, ut, dut, v) = dut * v

  t_Ω = FETerm(res, jac, jac_t, trian, quad)
  op = TransientFEOperator(U, V0, t_Ω)

  U0 = U(0.0)
  uh0 = interpolate_everywhere(U0, u(0.0))
  u0 = get_free_values(uh0)

  return op, u0

end

# Wrappers for residual and jacobian for DiffEqs
function diffeq_wrappers(op)

  ode_op = get_algebraic_operator(op)
  ode_cache = allocate_cache(ode_op)

  function residual(res, du, u, p, t)
    # TO DO (minor): Improve update_cache! st do nothing if same time t as in the cache
    # now it would be done twice (residual and jacobian)
    ode_cache = GridapODEs.ODETools.update_cache!(ode_cache, ode_op, t)
    GridapODEs.ODETools.residual!(res, ode_op, t, u, du, ode_cache)
  end

  function jacobian(jac, du, u, p, gamma, t)
    ode_cache = GridapODEs.ODETools.update_cache!(ode_cache, ode_op, t)
    z = zero(eltype(jac))
    Gridap.Algebra.fill_entries!(jac, z)
    GridapODEs.ODETools.jacobian_and_jacobian_t!(jac, ode_op, t, u, du, gamma, ode_cache)
  end

  function mass(mass, du, u, p, t)
    ode_cache = GridapODEs.ODETools.update_cache!(ode_cache, ode_op, t)
    z = zero(eltype(mass))
    Gridap.Algebra.fill_entries!(mass, z)
    GridapODEs.ODETools.jacobian_t!(mass, ode_op, t, u, du, 1.0, ode_cache)
  end

  function stiffness(stif, du, u, p, t)
    ode_cache = GridapODEs.ODETools.update_cache!(ode_cache, ode_op, t)
    z = zero(eltype(stif))
    Gridap.Algebra.fill_entries!(stif, z)
    GridapODEs.ODETools.jacobian!(stif, ode_op, t, u, du, ode_cache)
  end

  return residual, jacobian, mass, stiffness

end

# Allocate Jacobian
function allocate_jac(op,u0)
  ode_op = get_algebraic_operator(op)
  ode_cache = allocate_cache(ode_op) # Not acceptable in terms of performance
  return GridapODEs.ODETools.allocate_jacobian(ode_op, u0, ode_cache)
end

end #module


module DiffEqsWrapperTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator
using GridapODEs
using DifferentialEquations
using Gridap.Algebra: NewtonRaphsonSolver
using Base.Iterators

# Solving the heat equation using GridapODEs and DiffEqs
using Sundials
tspan = (0.0, 1.0)

u(x, t) = t
u(t) = x -> u(x, t)

# ISSUE 1: When I choose n > 2, even though the problem that we will solve is
# linear, the Sundials solvers seems to have convergence issues in the nonlinear
# solver (?). Ut returns errors
# [IDAS ERROR]  IDACalcIC Newton/Linesearch algorithm failed to converge.
# ISSUE 2: When I pass `jac_prototype` the code gets stuck
S = Main.GridapFESolver

n = 3 # cells per dim (2D)
op, u0 = S.fe_problem(u,n)

# Some checks
residual, jacobian, mass, stiffness = S.diffeq_wrappers(op)
J = S.allocate_jac(op,u0)
r = copy(u0)
θ = 1.0; t0 = 0.0; tF = 1.0; dt = 0.1; tθ = 1.0; dtθ = dt*θ
residual(r, u0, u0, nothing, tθ)
jacobian(J, u0, u0, nothing, (1 / dtθ), tθ)

K = S.allocate_jac(op,u0)
M = S.allocate_jac(op,u0)
stiffness(K, u0, u0, nothing, tθ)
mass(M, u0, u0, nothing, tθ)
# Here you have the mass matrix M

@test (1/dtθ)*M+K ≈ J

# To explore the Sundials solver options, e.g., BE with fixed time step dtd
f_iip = DAEFunction{true}(residual; jac = jacobian)#, jac_prototype=J)
# jac_prototype is the way to pass my pre-allocated jacobian matrix
prob_iip = DAEProblem{true}(f_iip, u0, u0, tspan, differential_vars = [true])
# When I pass `jac_prototype` the code get stuck here:
sol_iip = Sundials.solve(prob_iip, IDA(), reltol = 1e-8, abstol = 1e-8)
@show sol_iip.u

# or iterator version
integ = init(prob_iip, IDA(), reltol = 1e-8, abstol = 1e-8)
# step!(integ)

# Show using integrators as iterators
for i in take(integ, 100)
  @show integ.u
end

end # module

# FUTURE WORK: Check other options, not only Sundials

# ISSUE: Future work, add own (non)linear solver.
# Let us assume that we just want to consider a fixed point algorithm
# and we consider an implicit time integration of a nonlinear PDE.
# Our solvers are efficient since they re-use cache among
# nonlinear steps (e.g., symbolic factorization, arrays for numerical factorization)
# and for the transient case in our library, we can also reuse all this
# between time steps. Could we attain something like this using
# DifferentialEquations/Sundials?
# @ChrisRackauckas suggests to take a look at:
# https://docs.sciml.ai/latest/tutorials/advanced_ode_example/ shows swapping out linear solvers.
# https://docs.sciml.ai/latest/features/linear_nonlinear/ is all of the extra details.

# Try to pass solver too
# ls = LUSolver()
# ls_cache = nothing
# x = copy(u0)
# solve!(x,J,r,ls_cache)
#
