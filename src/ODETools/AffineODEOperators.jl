#
# """
# It represents an _affine_ operator in an implicit ODE, i.e., an ODE operator of
# the form A(t,u,∂tu) = M(t)∂tu + K(t)u + f(t)  where the
# implicit PDE reads A(t,u,∂tu) = 0, when ∂tu is the time derivative of u.
# The trait {C} stands for `is_constant`, if `true`, there is no dependency
# of M and K on time t.
# """
# abstract type AffineODEOperator <: ODEOperator end
#
# abstract type ConstantODEOperator <: AffineODEOperator end
#
# santiagobadia : These methods in a second step
# get_matrix(op::AffineODEOperator) = @abstractmethod
# get_matrix_t(op::AffineODEOperator) = @abstractmethod
# vector!(b,op::AffineODEOperator,t::Real,ode_cache) = @abstractmethod

#
# """
# It provides f(t) for a given t
# """
# function vector!(r::AbstractVector,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector,ode_cache)
#   @abstractmethod
# end
#
# """
# """
# function allocate_vector(op::ODEOperator,u::AbstractVector,ode_cache)
#   @abstractmethod
# end
#
# """
# It adds K(t) for a given (t) to a given matrix
# """
# function matrix!(J::AbstractMatrix,op::ODEOperator,t::Real,ode_cache)
#   @abstractmethod
# end
#
# """
# It adds γ*M(t) to a given matrix, where γ is as in the `ODEOperator`
# """
# function matrix_t!(J::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector,dut_u::Real,ode_cache)
#   @abstractmethod
# end
#
# """
# """
# function allocate_matrix(op::ODEOperator,u::AbstractVector,ode_cache)
#   @abstractmethod
# end
#
# """
# """
# function matrix_and_vector(op::ConstantAffineODEOperator,ode_cache::nothing)
#   A = allocate_matrix(op,)
#   A, b = cache
# end
#
# function matrix_and_vector(A,b,op::ConstantAffineODEOperator,)
#   A, b = cache
# end
