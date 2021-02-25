"""
A wrapper of `TransientFEOperator` that transforms it to `ODEOperator`, i.e.,
takes A(t,uh,∂tuh,vh) and returns A(t,uF,∂tuF) where uF and ∂tuF represent the
free values of the `EvaluationFunction` uh and ∂tuh.
"""
struct ODEOpFromFEOp{C} <: ODEOperator{C}
  feop::TransientFEOperator{C}
end

get_order(op::ODEOpFromFEOp) = get_order(op.feop)

function allocate_cache(op::ODEOpFromFEOp)
  Ut = get_trial(op.feop)
  U = allocate_trial_space(Ut)
  Uts = (Ut,)
  Us = (U,)
  for i in 1:get_order(op)
    Uts = (Uts...,∂t(Uts[i]))
    Us = (Us...,allocate_trial_space(Uts[i+1]))
  end
  fecache = allocate_cache(op.feop)
  ode_cache = (Us,Uts,fecache)
  ode_cache
end

function allocate_cache(op::ODEOpFromFEOp,v::AbstractVector,a::AbstractVector)
  ode_cache = allocate_cache(op)
  (v,a, ode_cache)
end

function update_cache!(ode_cache,op::ODEOpFromFEOp,t::Real)
  _Us,Uts,fecache = ode_cache
  Us = ()
  for i in 1:get_order(op)+1
    Us = (Us...,evaluate!(_Us[i],Uts[i],t))
  end
  fecache = update_cache!(fecache,op.feop,t)
  (Us,Uts,fecache)
end

function allocate_residual(op::ODEOpFromFEOp,uhF::AbstractVector,ode_cache)
  Us,Uts,fecache = ode_cache
  uh = EvaluationFunction(Us[1],uhF)
  allocate_residual(op.feop,uh,fecache)
end

function allocate_jacobian(op::ODEOpFromFEOp,uhF::AbstractVector,ode_cache)
  Us,Uts,fecache = ode_cache
  uh = EvaluationFunction(Us[1],uhF)
  allocate_jacobian(op.feop,uh,fecache)
end

function residual!(
  b::AbstractVector,
  op::ODEOpFromFEOp,
  t::Real,
  xhF::Tuple{Vararg{AbstractVector}},
  ode_cache)
  Xh, = ode_cache
  xh = ()
  for i in 1:get_order(op)+1
    xh = (xh...,EvaluationFunction(Xh[i],xhF[i]))
  end
  residual!(b,op.feop,t,xh,ode_cache)
end

function jacobian!(
  A::AbstractMatrix,
  op::ODEOpFromFEOp,
  t::Real,
  xhF::Tuple{Vararg{AbstractVector}},
  i::Integer,
  γᵢ::Real,
  ode_cache)
  Xh, = ode_cache
  xh = ()
  for i in 1:get_order(op)+1
    xh = (xh...,EvaluationFunction(Xh[i],xhF[i]))
  end
  jacobian!(A,op.feop,t,xh,i,γᵢ,ode_cache)
end

function jacobians!(
  J::AbstractMatrix,
  op::ODEOpFromFEOp,
  t::Real,
  xhF::Tuple{Vararg{AbstractVector}},
  γ::Tuple{Vararg{Real}},
  ode_cache)
  Xh, = ode_cache
  xh = ()
  for i in 1:get_order(op)+1
    xh = (xh...,EvaluationFunction(Xh[i],xhF[i]))
  end
  jacobians!(J,op.feop,t,xh,γ,ode_cache)
end
