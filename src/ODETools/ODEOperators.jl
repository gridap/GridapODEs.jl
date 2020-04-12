"""
It represents the operator in an implicit ODE, i.e., A(t,u,∂tu) where the
implicit PDE reads A(t,u,∂tu) = 0, when ∂tu is the time derivative of u
"""
abstract type ODEOperator <: GridapType end

# @santiagobadia : I would consider in a future a more general case, in which
# the implicit ODE has an arbitrary order, i.e., A(t,u,u_t, u_nt) = 0.
# get_order(op::ODEOperator) = @notimplemented
# @santiagobadia :
# We probably want to consider second order time derivatives too, i.e.,
# A(t,u,u_t,u_tt) = 0 (wave propagation, elastodynamics) or even more...

"""
It provides A(t,u,∂tu) for a given (t,u,∂tu)
"""
function residual!(r::AbstractVector,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector,op_state)
  @abstractmethod
end

"""
"""
function allocate_residual(op::ODEOperator,u::AbstractVector,state)
  @abstractmethod
end

"""
It adds [∂A/∂u](t,u,∂tu) for a given (t,u,∂tu) to a given matrix J
"""
function jacobian!(J::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector,op_state)
  @abstractmethod
  # Add values to J
end

"""
Add the contribution γ*[∂A/∂(∂tu)](u,∂tu) to the Jacobian matrix, where γ
is a scaling coefficient provided by the `ODESolver`, e.g., 1/Δt for Backward
Euler; It represents ∂(δt(u))/∂(u), in which δt(⋅) is the approximation of ∂t(⋅)
in the solver.
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

"""
Allocates the state data required by the `ODESolution` for a given `ODEOperator`
"""
allocate_state(op::ODEOperator) = @notimplemented

update_state!(state,op::ODEOperator,t::Real) = @notimplemented

"""
Tests the interface of `ODEOperator` specializations
"""
function test_ode_operator(op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  state = allocate_state(op)
  state = update_state!(state,op,0.0)
  r = allocate_residual(op,u,state)
  residual!(r,op,t,u,u_t,state)
  J = allocate_jacobian(op,u,state)
  jacobian!(J,op,t,u,u_t,state)
  jacobian_t!(J,op,t,u,u_t,1.0,state)
  true
end
