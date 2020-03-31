# EXAMPLE:

# Toy linear ODE with 2 DOFs
# u_1_t = a * u_1
# u_2_t = b * u_1 + c * u_2
struct ODEOperatorMock{T<:Real} <: ODEOperator
  a::T
  b::T
  c::T
end

function residual!(r::AbstractVector,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  r .= 0
  r[1] = - u_t[1] + op.a * u[1]
  r[2] = - u_t[2] + op.b * u[1] + op.c * u[2]
  r
end

function allocate_residual(op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  zeros(2)
end

function jacobian_unk!(j_u::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  j_u .= 0
  j_u[1,1] = op.a
  j_u[2,1] = op.b
  j_u[2,2] = op.c
  j_u
end

function jacobian_unk_t!(j_u_t::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  j_u_t .= 0
  j_u_t[1,1] = -1
  j_u_t[2,2] = -1
  j_u_t
end

function allocate_jacobian(op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector)
  zeros(2,2)
end
