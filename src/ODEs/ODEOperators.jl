# @santiagobadia Here we are defining A(t,u,u_t) in A(t,u,u_t) = 0.
abstract type ODEOperator <: GridapType end
# @santiagobadia: There is quite a lot of noise about what an ODE is, i.e.,
# whether it requires a non-singular Jacobian wrt u_t (DAE in this case) or not.
# I decided to use the term TimeStepper to circumvent this issue. In some places
# they use the term implicit ODEs to what I propose to use.

# For me the ODEOperator should be quite simple. In the Gridap FEM context,
# it will take a residual r(t,u,u_t,v), jacobian_unk(t,u,u_t,v) and
# jacobian_unk_t(t,u,u_t,v). In any case, there are some things to be
# considered when designing this object. TrialSpaces can also depend on time,
# due to time dependent boundary conditions. Also, new transient functions
# will be needed, e.g., u_D(x,t). We must think how to implement all this.
# We should implement a FEODEOperator (or similar name) for this

# get_order(op::ODEOperator) = @notimplemented
# @santiagobadia :
# We probably want to consider second order time derivatives too, i.e.,
# A(t,u,u_t,u_tt) = 0 (wave propagation, elastodynamics) or even more...

# This ones represents the value `A(t,u,v)`
function residual!(r::AbstractVector,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  @abstractmethod
end
# @santiagobadia : Agreed


function allocate_residual(op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  @abstractmethod
end
# @santiagobadia : Agreed

# This one represents `[∂A/∂u](t,u,v)`
function jacobian_unk!(j_u::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  @abstractmethod
end
# @santiagobadia : Agreed
# @santiagobadia : I think unk and unk_t or u u_t is more revealing than u v
# Not exposed to user anyway

# This one represents `[∂A/∂v](t,u,v)`
function jacobian_unk_t!(j_u_t::AbstractMatrix,op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  @abstractmethod
end
# @santiagobadia : Agreed

function allocate_jacobian(op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  @abstractmethod
end
# @santiagobadia : Agreed
