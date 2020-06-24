"""
Trait for `ODEOperator` that tells us whether the operator depends on the solution
(including its time derivatives), it is an affine operator that depends on time
or it is a constant operator (affine and time-indepedendent)
"""
abstract type OperatorType end
struct Nonlinear <: OperatorType end
struct Affine  <: OperatorType end
struct Constant  <: OperatorType end

"""
It represents the operator in an implicit ODE, i.e., A(t,u,∂tu) where the
implicit PDE reads A(t,u,∂tu) = 0, when ∂tu is the time derivative of u.
The trait `{C}` determines whether the operator is fully nonlinear, affine
or constant in time.
"""
abstract type ODEOperator{C<:OperatorType} <: GridapType end

"""
It represents an _affine_ operator in an implicit ODE, i.e., an ODE operator of
the form A(t,u,∂tu) = M(t)∂tu + K(t)u + f(t)
"""
const AffineODEOperator = ODEOperator{Affine}
# abstract type AffineODEOperator <: ODEOperator end

"""
It represents a constant operator in an implicit ODE, i.e., an ODE operator of
the form A(t,u,∂tu) = M∂tu + Ku + f
"""
const ConstantODEOperator = ODEOperator{Constant}
# abstract type ConstantODEOperator <: AffineODEOperator end

# @santiagobadia : I would consider in a future a more general case, in which
# the implicit ODE has an arbitrary order, i.e., A(t,u,u_t, u_nt) = 0.
# get_order(op::ODEOperator) = @notimplemented
# @santiagobadia :
# We probably want to consider second order time derivatives too, i.e.,
# A(t,u,u_t,u_tt) = 0 (wave propagation, elastodynamics) or even more...
"""
Returns the `OperatorType`, i.e., nonlinear, affine, or constant in time
"""
OperatorType(::ODEOperator{C}) where C = C
# OperatorType(::Type{<:ODEOperator{C}}) where C = C

"""
It provides A(t,u,∂tu) for a given (t,u,∂tu)
"""
function residual!(r::AbstractVector,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector,ode_cache)
  @abstractmethod
end

"""
"""
function allocate_residual(op::ODEOperator,u::AbstractVector,ode_cache)
  @abstractmethod
end

"""
It adds [∂A/∂u](t,u,∂tu) for a given (t,u,∂tu) to a given matrix J
"""
function jacobian!(J::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector,ode_cache)
  @abstractmethod
  # Add values to J
end

"""
Add the contribution γ*[∂A/∂(∂tu)](u,∂tu) to the Jacobian matrix, where γ
is a scaling coefficient provided by the `ODESolver`, e.g., 1/Δt for Backward
Euler; It represents ∂(δt(u))/∂(u), in which δt(⋅) is the approximation of ∂t(⋅)
in the solver.
"""
function jacobian_t!(J::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector,dut_u::Real,ode_cache)
  @abstractmethod
  # Add values to J
end

"""
"""
function allocate_jacobian(op::ODEOperator,u::AbstractVector,ode_cache)
  @abstractmethod
end

"""
Allocates the cache data required by the `ODESolution` for a given `ODEOperator`
"""
allocate_cache(op::ODEOperator) = @abstractmethod

#@fverdugo to be used as `cache = update_cache!(cache,op,t)`
update_cache!(cache,op::ODEOperator,t::Real) = @abstractmethod

"""
Tests the interface of `ODEOperator` specializations
"""
function test_ode_operator(op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  cache = allocate_cache(op)
  cache = update_cache!(cache,op,0.0)
  r = allocate_residual(op,u,cache)
  residual!(r,op,t,u,u_t,cache)
  J = allocate_jacobian(op,u,cache)
  jacobian!(J,op,t,u,u_t,cache)
  jacobian_t!(J,op,t,u,u_t,1.0,cache)
  true
end
