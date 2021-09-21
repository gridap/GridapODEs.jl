# GridapODEs

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/GridapODEs.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/GridapODEs.jl/dev)
[![Build Status](https://github.com/gridap/GridapODEs.jl/workflows/CI/badge.svg?branch=master)](https://github.com/gridap/GridapODEs.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/gridap/GridapODEs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/GridapODEs.jl)
[![Coveralls](https://coveralls.io/repos/github/gridap/GridapODEs.jl/badge.svg?branch=master)](https://coveralls.io/github/gridap/GridapODEs.jl?branch=master)
[![DOI](https://zenodo.org/badge/250735390.svg)](https://zenodo.org/badge/latestdoi/250735390)

This package provides time integration tools for `Gridap`. As an example, the following code solves the heat equation.

```julia
using Gridap
using ForwardDiff
using GridapODEs.ODETools
using GridapODEs.TransientFETools

import Gridap: ∇
import GridapODEs.TransientFETools: ∂t

θ = 0.5

u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(t) = x -> u(x,t)

f(t) = x -> ∂t(u)(x,t)-Δ(u(t))(x) # or ∂t(u)(t)(x)-Δ(u(t))(x)

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

order = 2

reffe = ReferenceFE(lagrangian,Float64,order)
V0 = FESpace(
  model,
  reffe,
  conformity=:H1, 
  dirichlet_tags="boundary"
)

U = TransientTrialFESpace(V0,u)

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)

a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
b(v,t) = ∫( v*f(t) )dΩ
m(u,v) = ∫( v*u )dΩ

res(t,u,ut,v) = a(u,v) + m(ut,v) - b(v,t)
jac(t,u,ut,du,v) = a(du,v)
jac_t(t,u,ut,dut,v) = m(dut,v)

op = TransientFEOperator(res,jac,jac_t,U,V0)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)

ls = LUSolver()
odes = ThetaMethod(ls,dt,θ)
solver = TransientFESolver(odes)
sol_t = solve(solver,op,uh0,t0,tF)

for (uh_tn, tn) in sol_t
  # Here we have the solution uh_tn at tn
end
```
