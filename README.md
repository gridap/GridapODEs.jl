# GridapODEs

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/GridapODEs.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/GridapODEs.jl/dev)
[![Build Status](https://travis-ci.com/gridap/GridapODEs.jl.svg?branch=master)](https://travis-ci.com/gridap/GridapODEs.jl)
[![Codecov](https://codecov.io/gh/gridap/GridapODEs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/GridapODEs.jl)
[![Coveralls](https://coveralls.io/repos/github/gridap/GridapODEs.jl/badge.svg?branch=master)](https://coveralls.io/github/gridap/GridapODEs.jl?branch=master)

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
u(t::Real) = x -> u(x,t)
∂tu(t) = x -> ForwardDiff.derivative(v(x),t)
∂tu(x,t) = ∂tu(t)(x)

f(t) = x -> 1.0

domain = (0,1,0,1)
partition = (4,4)
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
solver = TransientFESolver(odes)
sol_t = solve(solver,op,uh0,t0,tF)

for (uh_tn, tn) in sol_t
  # Here we have the solution uh_tn at tn
end
```
