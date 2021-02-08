# Toy linear ODE with 2 DOFs
# u_1_t - a * u_1 = 0
# u_2_t - b * u_1 - c * u_2 = 0

# Toy 2nd order ODE with 2 DOFs
# u_1_tt + b * u_1_t - a * u_1 = 0
# u_2_tt + a * u_1_t - b * u_1 - c * u_2 = 0

import GridapODEs.ODETools: ODEOperator
import GridapODEs.ODETools: AffineODEOperator
import GridapODEs.ODETools: ConstantODEOperator
import GridapODEs.ODETools: allocate_cache
import GridapODEs.ODETools: update_cache!
import GridapODEs.ODETools: allocate_residual
import GridapODEs.ODETools: jacobian!
import GridapODEs.ODETools: jacobian_t!
import GridapODEs.ODETools: jacobian_and_jacobian_t!
import GridapODEs.ODETools: allocate_jacobian
import GridapODEs.ODETools: residual!

struct ODEOperatorMock{T<:Real,C} <: ODEOperator{C}
  a::T
  b::T
  c::T
end

function residual!(r::AbstractVector,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector,ode_cache)
  r .= 0
  r[1] = u_t[1] - op.a * u[1]
  r[2] = u_t[2] - op.b * u[1] - op.c * u[2]
  r
end

function residual!(r::AbstractVector,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector,u_tt::AbstractVector,ode_cache)
  r .= 0
  r[1] = u_tt[1] + op.b * u_t[1] - op.a * u[1]
  r[2] = u_tt[2] + op.a * u_t[1]- op.b * u[1] - op.c * u[2]
  r
end

function allocate_residual(op::ODEOperatorMock,u::AbstractVector,cache)
  zeros(2)
end

function jacobian!(J::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector,ode_cache)
  J[1,1] += -op.a
  J[2,1] += -op.b
  J[2,2] += -op.c
  J
end

function jacobian!(J::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector,u_tt::AbstractVector,ode_cache)
  J[1,1] += -op.a
  J[2,1] += -op.b
  J[2,2] += -op.c
  J
end

function jacobian_t!(J::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector,du_t_u::Real,ode_cache)
  J[1,1] += 1.0*du_t_u
  J[2,2] += 1.0*du_t_u
  J
end

function jacobian_t!(J::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector,u_tt::AbstractVector,du_t_u::Real,ode_cache)
  J[1,1] += op.b*du_t_u
  J[2,1] += op.a*du_t_u
  J
end

function jacobian_tt!(J::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector,u_tt::AbstractVector,du_tt_u::Real,ode_cache)
  J[1,1] += 1.0*du_tt_u
  J[2,2] += 1.0*du_tt_u
  J
end

function jacobian_and_jacobian_t!(J::AbstractMatrix,op::ODEOperatorMock,t::Real,u::AbstractVector,u_t::AbstractVector,du_t_u::Real,ode_cache)
  jacobian!(J,op,t,u,u_t,ode_cache)
  jacobian_t!(J,op,t,u,u_t,du_t_u,ode_cache)
  J
end

function jacobian_and_jacobian_t!(
  J::AbstractMatrix,
  op::ODEOperatorMock,
  t::Real,
  u::AbstractVector,
  u_t::AbstractVector,
  u_tt::AbstractVector,
  du_t_u::Real,
  du_tt_u::Real,
  ode_cache)
  jacobian!(J,op,t,u,u_t,u_tt,ode_cache)
  jacobian_t!(J,op,t,u,u_t,u_tt,du_t_u,ode_cache)
  jacobian_tt!(J,op,t,u,u_t,u_tt,du_tt_u,ode_cache)
  J
end



function allocate_jacobian(op::ODEOperatorMock,u::AbstractVector,cache)
  zeros(2,2)
end

allocate_cache(op::ODEOperatorMock) = nothing
allocate_cache(op::ODEOperatorMock,v::AbstractVector,a::AbstractVector) = (v,a),nothing
update_cache!(cache,op::ODEOperatorMock,t::Real) = cache
