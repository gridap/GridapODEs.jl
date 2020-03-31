module GridapTimeStepper
```
sdsdasda
dasdasdas
```
using Gridap
# For the moment I am going to use this file for brainstorming..

abstract type TimeStepper end

# It is an iterator, thus we need to implement the minimum
# Iteration interface for it to work
# see https://docs.julialang.org/en/v1/manual/interfaces/

iterate(::TimeStepper) = @notimplemented

iterate(::TimeStepper, state) = @notimplemented

# Other optional methods that can be implemented by
# concrete types

# IteratorSize(::TimeStepper) = HasSize()
# Or HasLength(), HasShape{N}(), IsInfinite(), or SizeUnknown()

# IteratorEltype(::TimeStepper) = HasEltype()
# Replace with HasEltype() for concrete TimeStepper

# eltype(::TimeStepper)	= Any

# length(::TimeStepper) = undef

# size(::TimeStepper,dim) = undef

```
A system of ODE/DAE equations can be written as

A(t,u,u_t) = 0.

Given a point p = (t_,u_,∂t_u_), we need to provide methods that evaluate

A(p), ∂A/∂u(p), ∂A/∂u_t(p)

The `TimeStepper` will have internal information that will allow it to do the
time integration. E.g., for a θ-method, given u0

uθ = θuF + (1-θ)u0, or uF = (1/θ)*uθ - ((1-θ)/θ)*u0, θ > 0,
u_t = (uF-u0)/δt = (uθ-u0)/(θδt),
∂u_t/∂uθ = 1/(δt*θ),

and the nonlinear problem to be solved has the Jacobian

(∂A/∂u*(1/θ) + ∂A/∂u_t(p))(tθ_,uθ_,u_t_) δuθ = -A(tθ_,uθ_,u_t_).

Thus, we can compute the map TSθ(u0) = uF.

So, I guess we need a @law that will accept time too, e.g., @transientlaw

Anything else?

The only question is what we are going to provide to the `TimeStepper` in order
to proceed. I guess that we should create a new `FETerm`, e.g.,
`TransientFETerm` that requires the residual and the Jacobians with respect to
both `u` and `u_t` and a `TransientFEOperator` for this transient term (not sure
about the last one).

On the other hand, we should provide `TimeStepper` constructor with a
(non)linear solver in general. That seems easy the way we have structured things.

```


end # module

# Like we have done in the spatial FE discretization, I would distinguish between the
# algebraic and FE counterparts and between the problem to be solved and the method to solve the problem.
#
# Disclaimer: The names of functions and types I use here are just an initial guess. They can be improved a lot.
#
# For the moment, we can focus in the algebraic part.
# From this, designing the FE counterpart will be easy.
#
# First ingredient: object describing the algebraic problem to solve.
# We can follow the notation A(t,u,u_t)=0 you have proposed to represent the ODE.
# For simplicity, I use here `A(t,u,v)=0`.
# The first we need is an abstract type representing the operator A:
#
