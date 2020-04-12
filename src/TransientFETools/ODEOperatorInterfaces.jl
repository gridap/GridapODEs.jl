# A concrete version of ODEOperator from a FEOperator
# struct ODEOperatorFromFEOperator <: ODEOperator
#   op::FEOperator
# end

struct ODEOpFromFEOp <: ODEOperator
  feop::TransientFEOperator
end

function allocate_cache(op::ODEOpFromFEOp)
  allocate_cache(op.feop)
end

function update_cache!(cache,op::ODEOpFromFEOp,t::Real)
  update_cache!(cache,op.feop,t)
end

function allocate_residual(op::ODEOpFromFEOp,uhF::AbstractVector,ode_cache)
  Uh = ode_cache.Uh; Uht = ode_cache.Uht
  uh = FEFunction(Uh,uhF)
  allocate_residual(op.feop,uh)#,ode_cache)
end

function allocate_jacobian(op::ODEOpFromFEOp,uhF::AbstractVector,ode_cache)
  Uh = ode_cache.Uh; Uht = ode_cache.Uht
  uh = FEFunction(Uh,uhF)
  allocate_jacobian(op.feop,uh,ode_cache)
end

function residual!(b::AbstractVector,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,ode_cache)
  Uh = ode_cache.Uh; Uht = ode_cache.Uht
  uh = FEFunction(Uh,uhF)
  uht = FEFunction(Uht,uhtF)
  residual!(b,op.feop,t,uh,uht,ode_cache)
end

function jacobian!(A::AbstractMatrix,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,ode_cache)
  Uh = ode_cache.Uh; Uht = ode_cache.Uht
  uh = FEFunction(Uh,uhF)
  uht = FEFunction(Uht,uhtF)
  jacobian!(A,op.feop,t,uh,uht,ode_cache)
end

function jacobian_t!(J::AbstractMatrix,op::ODEOpFromFEOp,t::Real,uhF::AbstractVector,uhtF::AbstractVector,dut_u::Real,ode_cache)
  Uh = ode_cache.Uh; Uht = ode_cache.Uht
  uh = FEFunction(Uh,uhF)
  uht = FEFunction(Uht,uhtF)
  jacobian_t!(J,op.feop,t,uh,uht,dut_u,ode_cache)
end
