# A concrete version of ODEOperator from a FEOperator
# struct ODEOperatorFromFEOperator <: ODEOperator
#   op::FEOperator
# end
"""
A wrapper of `TransientFEOperator` that transforms it to `ODEOperator`, i.e.,
takes A(t,uh,∂tuh,vh) and returns A(t,uF,∂tuF) where uF and ∂tuF represent the
free values of the `FEFunction` uh and ∂tuh.
"""
struct ODEOpFromFEOp <: ODEOperator
  feop::TransientFEOperator
end

function allocate_cache(op::ODEOpFromFEOp)
  Ut = get_trial(op.feop)
  Vt = ∂t(Ut)
  U = allocate_trial_space(Ut)
  V = allocate_trial_space(Vt)
  fecache = allocate_cache(op.feop)
  ode_cache = (U,V,Ut,Vt,fecache)
  ode_cache
end

function update_cache!(ode_cache,op::ODEOpFromFEOp,t::Real)
  U,V,Ut,Vt,fecache = ode_cache
  U = evaluate!(U,Ut,t)
  V = evaluate!(V,Vt,t)
  fecache = update_cache!(fecache,op.feop,t)
  (U,V,Ut,Vt,fecache)
end

function allocate_residual(op::ODEOpFromFEOp,uhF::AbstractVector,ode_cache)
  U,V,Ut,Vt,fecache = ode_cache
  uh = FEFunction(Uh,uhF)
  allocate_residual(op.feop,uh,fecache)
end

function allocate_jacobian(op::ODEOpFromFEOp,uhF::AbstractVector,ode_cache)
  U,V,Ut,Vt,fecache = ode_cache
  uh = FEFunction(Uh,uhF)
  allocate_jacobian(op.feop,uh,fecache)
end

function residual!(b::AbstractVector,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,ode_cache)
  Uh,Uht, = ode_cache
  uh = FEFunction(Uh,uhF)
  uht = FEFunction(Uht,uhtF)
  residual!(b,op.feop,t,uh,uht,ode_cache)
end

function jacobian!(A::AbstractMatrix,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,ode_cache)
  Uh,Uht, = ode_cache
  uh = FEFunction(Uh,uhF)
  uht = FEFunction(Uht,uhtF)
  jacobian!(A,op.feop,t,uh,uht,ode_cache)
end

function jacobian_t!(J::AbstractMatrix,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,dut_u::Real,ode_cache)
  Uh,Uht, = ode_cache
  uh = FEFunction(Uh,uhF)
  uht = FEFunction(Uht,uhtF)
  jacobian_t!(J,op.feop,t,uh,uht,dut_u,ode_cache)
end
