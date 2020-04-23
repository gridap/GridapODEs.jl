using DifferentialEquations

f(u,p,t) = 1.01*u
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, ImplicitEuler(), reltol=1e-8, abstol=1e-8)

using Plots
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")


# DAE problem
u0 = zeros(1)
u0[1] = 1/2

function residual(res,du,u,p,t)
  res[1] = - 1.01u[1] + du[1]
end

function jacobian(jac,du,u,p,gamma,t)
  jac[1,1] = gamma - 1.01
  nothing
end

function residual_ode(du,u,p,t)
  du[1] = 1.01u[1]
end

function update_func(mass,u,p,t)
  mass = ones(1,1)
end

dependent_M = DiffEqArrayOperator(ones(1,1),update_func=update_func)

# Using https://docs.sciml.ai/v2.1.0/features/performance_overloads.html#Declaring-Explicit-Jacobians-1
function residual_ode(::Type{Val{:jac}},t,u,J)
  jac[1,1] = - 1.01
  nothing
end

f_iip = ODEFunction{true}(residual_ode, mass_matrix=dependent_M)
prob_iip = ODEProblem{true}(f_iip,u0,tspan)
sol_iip = solve(prob_iip, ImplicitEuler(), reltol=1e-8, abstol=1e-8)

# f_iip = ODEFunction{true}(residual;jac=jacobian)
# Look here for more options
#https://docs.sciml.ai/latest/features/performance_overloads/#
# How can we provide a solver for W = M - gamma*J or W_t = gamma*W


# Why is it looking for the OOP version of the Jacobian?

Thanks for the quick reply, I will consider this option.

In any case, following

https://docs.sciml.ai/latest/features/performance_overloads/#

I am creating a DAE problem

```julia
tspan = (0.0,1.0)
u0 = zeros(1)
u0[1] = 1/2

function residual(res,du,u,p,t)
  res[1] = - 1.01u[1] + du[1]
end

function jacobian(jac,du,u,p,gamma,t)
  jac[1,1] = gamma - 1.01
  nothing
end

f_iip = DAEFunction{true}(residual;jac=jacobian)
prob_iip = DAEProblem{true}(f_iip,u0,u0,tspan)
sol_iip = solve(prob_iip, DImplicitEuler(), reltol=1e-8, abstol=1e-8)
```
but it returns an error,
```
type ODEIntegrator has no field duprev
```
which is weird, since I would expect a `DAEIntegrator` being created.

It would be great if you could tell me what I am doing wrong.

Thanks

# f_iip = ODEFunction{true}(residual)#,mass_matrix=M)
# prob_iip = ODEProblem{true}(f_iip,u0,tspan)
# sol_iip = solve(prob_iip, Rodas5(), reltol=1e-8, abstol=1e-8)

# Juno.@enter solve(prob_iip, DImplicitEuler(), reltol=1e-8, abstol=1e-8)
# I understand the problem, in place overwrites the du array, BUT not in terms of residual









I have the following questions.

I think that even if my "mass" matrix is non-singular, I would like to use the
ADE interface you provide, since the mass matrix can be nonlinear, and thus, I
need to compute its derivative with respect to `du` in the Jacobian. The DAE
interface in `DifferentialEquations.jl` seems to provide this option and
perfectly matches the interface in our PDE library `Gridap`.

I would like to know whether using such interface instead of a ODE interface
has an impact on performance (assuming we can solve the same problem in both
cases).

######

I am trying to combine the ODE solvers in `DifferentialEquations.jl` with the PDE solvers in the `Gridap.jl` library.

It is mandatory for our applications to use the inplace version of the `DifferentialEquations.jl` package. This way, we can pre-allocate the required arrays in our `Gridap` solvers.

On the other hand, even when considering ODE systems, we want to use an implicit ODE expression of the solver, i.e., `A(du,u,p,t)=0`, since our _mass_ matrices are not just identity matrices, and can be fairly complicated or even nonlinear.

Finally, what we would provide is a `residual(res,du,u,p,t)` and `jacobian(J,du,u,p,t)` functions.

As a result, I think that the way to go is to consider an IIP problem in `DifferentialEquations.jl`. It seems that only the `DAE` constructor allows the syntax before.

However, for implicit ODE systems I cannot use an `ODE` algorithm, e.g., Backward-Euler, even though it is acceptable for our problems when the _mass_ matrix is non-singular (thus not a `DAE`). And it seems we must provide an initial value for `du` that does not have mathematical sense in our problems.

Summarizing, is there a way to create an `ODE` solver that takes `residual!(res,du,u,p,t)` and `jacobian!(J,du,u,p,t)` inplace functions, and allows me to use a `ODE` algorithm?

I have written the following silly linear example with what I would like to have

```julia
tspan = (0.0,1.0)
u0 = zeros(1)
u0[1] = 1/2

function residual(res,du,u,p,t)
  res[1] = - 1.01u[1] + du[1]
end

function jacobian(jac,du,u,p,gamma,t)
  jac[1,1] = gamma - 1.01
  nothing
end

f_iip = ODEFunction{true}(residual;jac=jacobian)
prob_iip = ODEProblem{true}(f_iip,u0,tspan)
sol_iip = solve(prob_iip, ImplicitEuler(), reltol=1e-8, abstol=1e-8)
```

but it does not work, because it seems the library is looking for an `ODE`-style inplace version of `residual`, i.e., `residual(du,u,p,t)`, based on the error:

```
MethodError: no method matching residual(::Array{Float64,1}, ::Array{Float64,1}, ::DiffEqBase.NullParameters, ::Float64)
Closest candidates are:
  residual(::Any, ::Any, ::Any, ::Any, !Matched::Any) at /.../DifEqWrappers.jl:20
...
```

In fact, even without considering the Jacobian, the same problem appears with the residual when doing this:
```julia
prob_iip = ODEProblem{false}(residual,u0,tspan)
sol_iip = solve(prob_iip, ImplicitEuler(), reltol=1e-8, abstol=1e-8)
```

Is it possible to do what I want in `DifferentialEquations.jl`? If Yes, What am I doing wrong? I could not find an example in the documentation covering this (probably my fault).

# On the other hand, is there a version of the solver that returns an iterator over the sequence of solutions in time?
#
# It would be a good feature, in complex PDEs one probably wants to consider the solution at each time step, perform
# space, compute a norm, etc, but not to store the solution at all time steps at
# the same time, since it is a lot of memory for a large PDE problem
