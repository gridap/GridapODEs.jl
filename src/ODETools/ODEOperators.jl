"""
It represents the operator in an implicit ODE, i.e., A(t,u,u_t) where the
implicit PDE reads A(t,u,u_t) = 0.
"""
abstract type ODEOperator <: GridapType end

# @santiagobadia : I would consider in a future a more general case, in which
# the implicit ODE has an arbitrary order, i.e., A(t,u,u_t, u_nt) = 0.
# get_order(op::ODEOperator) = @notimplemented
# @santiagobadia :
# We probably want to consider second order time derivatives too, i.e.,
# A(t,u,u_t,u_tt) = 0 (wave propagation, elastodynamics) or even more...

"""
It provides A(t,u,u_t) for a given (t,u,u_t)
"""
function residual!(r::AbstractVector,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector,op_state)
  @abstractmethod
end

function allocate_residual(op::ODEOperator,u::AbstractVector,state)
  @abstractmethod
end

"""
It adds [∂A/∂u](t,u,u_t) for a given (t,u,u_t) to a given matrix J
"""
function jacobian!(J::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector,op_state)
  @abstractmethod
  # Add values to J
end

"""
It adds [∂A/∂u_t](t,u,u_t) for a given (t,u,u_t) to a given matrix J
"""
function jacobian_t!(J::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector,dut_u::Real,op_state)
  @abstractmethod
  # Add values to J
end

"""
"""
function allocate_jacobian(op::ODEOperator,u::AbstractVector)
  @abstractmethod
end

allocate_state(op::ODEOperator) = @notimplemented

update_state!(state,op::ODEOperator,t::Real) = @notimplemented

"""
"""
function test_ode_operator(op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  state = allocate_state(op)
  state = update_state!(state,op,0.0)
  r = allocate_residual(op,u,state)
  residual!(r,op,t,u,u_t,state)
  J = allocate_jacobian(op,u,state)
  jacobian!(J,op,t,u,u_t,state)
  jacobian_t!(J,op,t,u,u_t,state)
  true
end
