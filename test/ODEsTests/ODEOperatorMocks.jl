# Toy linear ODE with 2 DOFs
# u_1_t - a * u_1 = 0
# u_2_t - b * u_1 - c * u_2 = 0

import GridapTimeStepper.ODETools: ODEOperator
import GridapTimeStepper.ODETools: allocate_residual
import GridapTimeStepper.ODETools: jacobian_unknown!
import GridapTimeStepper.ODETools: jacobian_unknown_t!
import GridapTimeStepper.ODETools: allocate_jacobian
import GridapTimeStepper.ODETools: residual!

struct ODEOperatorMock{T<:Real} <: ODEOperator
  a::T
  b::T
  c::T
end

function residual!(r::AbstractVector,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  r .= 0
  r[1] = u_t[1] - op.a * u[1]
  r[2] = u_t[2] - op.b * u[1] - op.c * u[2]
  r
end

function allocate_residual(op::ODEOperatorMock,u::AbstractVector,u_t::AbstractVector)
  zeros(2)
end

function jacobian_unknown!(J::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  J[1,1] += op.a
  J[2,1] += op.b
  J[2,2] += op.c
  J
end

function jacobian_unknown_t!(J::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  J[1,1] += 1.0
  J[2,2] += 1.0
  J
end

function allocate_jacobian(op::ODEOperatorMock,u::AbstractVector,u_t::AbstractVector)
  zeros(2,2)
end
