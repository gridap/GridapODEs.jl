module ODEOperatorsTests

using GridapTimeStepper.ODETools
using Test

import GridapTimeStepper.ODETools: test_ode_operator

include("ODEOperatorMocks.jl")

op = ODEOperatorMock(1.0,2.0,3.0)

u = ones(2)
u_t = ones(2)*2.0

@assert(length(u) == 2)
@assert(length(u_t) == 2)

r = allocate_residual(op,u,u_t)
@test r == zeros(2)

J = allocate_jacobian(op,u,u_t)
@test J == zeros(2,2)

t = 0.0
residual!(r,op,t,u,u_t)
_r = zeros(2)
_r[1] = u_t[1] - op.a * u[1]
_r[2] = u_t[2] - op.b * u[1] - op.c * u[2]
@test all(r .== _r)

J .= 0
jacobian_unknown!(J,op,t,u,u_t)
_J = zeros(2,2)
_J[1,1] = op.a
_J[2,1] = op.b
_J[2,2] = op.c
@test all(J .== _J)
_J
J
jacobian_unknown_t!(J,op,t,u,u_t)
_J[1,1] += 1.0
_J[2,2] += 1.0
@test all(J .== _J)
_J
J

@test test_ode_operator(op,t,u,u_t)

end #module
