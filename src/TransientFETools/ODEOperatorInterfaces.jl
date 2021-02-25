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
  Uts = (Ut)
  for i in 1:get_order(op)
    Uts = (Uts...,∂t(Uts[i]))
  end
  U = allocate_trial_space(Ut)
  Us = NTuple{get_order(op)+1,typeof(U)}
  for i in 1:get_order(op)
    Us[i+1] = allocate_trial_space(Uts[i+1])
  end
  fecache = allocate_cache(op.feop)
  ode_cache = (Us,Uts,fecache)
  ode_cache
end

function update_cache!(ode_cache,op::ODEOpFromFEOp,t::Real)
  Us,Uts,fecache = ode_cache
  for i in 1:get_order(op)+1
    Us[i] = evaluate!(Us[i],Uts[i],t)
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
  uh = EvaluationFunction(U,uhF)
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


# # """
# # A wrapper of `Transient2ndOrderFEOperator` that transforms it to `ODEOperator`, i.e.,
# # takes A(t,uh,∂tuh,∂ttuh,vh) and returns A(t,uF,∂tuF,∂ttuF) where uF, ∂tuF and ∂ttuF represent the
# # free values of the `EvaluationFunction` uh, ∂tuh and ∂ttuh.
# # """
# struct SecondOrderODEOpFromFEOp{C} <: SecondOrderODEOperator{C}
#   feop::TransientFEOperator{C}
# end

# function allocate_cache(op::SecondOrderODEOpFromFEOp,v::AbstractVector,a::AbstractVector)
#   Ut = get_trial(op.feop)
#   Vt = ∂t(Ut)
#   At = ∂tt(Ut)
#   U = allocate_trial_space(Ut)
#   V = allocate_trial_space(Vt)
#   A = allocate_trial_space(At)
#   fecache = allocate_cache(op.feop)
#   ode_cache = (U,V,A,Ut,Vt,At,fecache)
#   (v,a, ode_cache)
# end

# function update_cache!(ode_cache,op::SecondOrderODEOpFromFEOp,t::Real)
#   v,a, (U,V,A,Ut,Vt,At,fecache) = ode_cache
#   U = evaluate!(U,Ut,t)
#   V = evaluate!(V,Vt,t)
#   A = evaluate!(A,At,t)
#   fecache = update_cache!(fecache,op.feop,t)
#   (v,a,(U,V,A,Ut,Vt,At,fecache))
# end

# function allocate_residual(op::SecondOrderODEOpFromFEOp,uhF::AbstractVector,ode_cache)
#   U,V,A,Ut,Vt,At,fecache = ode_cache
#   uh = EvaluationFunction(U,uhF)
#   allocate_residual(op.feop,uh,fecache)
# end

# function allocate_jacobian(op::SecondOrderODEOpFromFEOp,uhF::AbstractVector,ode_cache)
#   U,V,A,Ut,Vt,At,fecache = ode_cache
#   uh = EvaluationFunction(U,uhF)
#   allocate_jacobian(op.feop,uh,fecache)
# end

# function residual!(
#   b::AbstractVector,
#   op::SecondOrderODEOpFromFEOp,
#   t::Real,
#   uhF::AbstractVector,
#   uhtF::AbstractVector,
#   uhttF::AbstractVector,
#   ode_cache)
#   Uh,Uht,Uhtt, = ode_cache
#   uh = EvaluationFunction(Uh,uhF)
#   uht = EvaluationFunction(Uht,uhtF)
#   uhtt = EvaluationFunction(Uhtt,uhttF)
#   residual!(b,op.feop,t,uh,uht,uhtt,ode_cache)
# end

# function jacobian!(
#   A::AbstractMatrix,
#   op::SecondOrderODEOpFromFEOp,
#   t::Real,
#   uhF::AbstractVector,
#   uhtF::AbstractVector,
#   uhttF::AbstractVector,
#   ode_cache)
#   Uh,Uht,Uhtt, = ode_cache
#   uh = EvaluationFunction(Uh,uhF)
#   uht = EvaluationFunction(Uht,uhtF)
#   uhtt = EvaluationFunction(Uhtt,uhttF)
#   jacobian!(A,op.feop,t,uh,uht,uhtt,ode_cache)
# end

# function jacobian_t!(
#   J::AbstractMatrix,
#   op::SecondOrderODEOpFromFEOp,
#   t::Real,
#   uhF::AbstractVector,
#   uhtF::AbstractVector,
#   uhttF::AbstractVector,
#   dut_u::Real,
#   ode_cache)
#   Uh,Uht,Uhtt, = ode_cache
#   uh = EvaluationFunction(Uh,uhF)
#   uht = EvaluationFunction(Uht,uhtF)
#   uhtt = EvaluationFunction(Uhtt,uhttF)
#   jacobian_t!(J,op.feop,t,uh,uht,uhtt,dut_u,ode_cache)
# end

# function jacobian_tt!(
#   J::AbstractMatrix,
#   op::SecondOrderODEOpFromFEOp,
#   t::Real,
#   uhF::AbstractVector,
#   uhtF::AbstractVector,
#   uhttF::AbstractVector,
#   dutt_u::Real,
#   ode_cache)
#   Uh,Uht,Uhtt, = ode_cache
#   uh = EvaluationFunction(Uh,uhF)
#   uht = EvaluationFunction(Uht,uhtF)
#   uhtt = EvaluationFunction(Uhtt,uhttF)
#   jacobian_t!(J,op.feop,t,uh,uht,uhtt,dutt_u,ode_cache)
# end

# function jacobian_and_jacobian_t!(
#   J::AbstractMatrix,
#   op::SecondOrderODEOpFromFEOp,
#   t::Real,
#   uhF::AbstractVector,
#   uhtF::AbstractVector,
#   uhttF::AbstractVector,
#   dut_u::Real,
#   dutt_u::Real,
#   ode_cache)
#   Uh,Uht,Uhtt, = ode_cache
#   uh = EvaluationFunction(Uh,uhF)
#   uht = EvaluationFunction(Uht,uhtF)
#   uhtt = EvaluationFunction(Uhtt,uhttF)
#   jacobian_and_jacobian_t!(J,op.feop,t,uh,uht,uhtt,dut_u,dutt_u,ode_cache)
# end
