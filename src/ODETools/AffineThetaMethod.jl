function solve_step!(uf::AbstractVector,
                     solver::ThetaMethod,
                     op::AffineODEOperator,
                     u0::AbstractVector,
                     t0::Real,
                     cache) # -> (uF,tF)

  dt = solver.dt
  solver.θ == 0.0 ? dtθ = dt : dtθ = dt*solver.θ
  tθ = t0+dtθ

  if cache === nothing
    ode_cache = allocate_cache(op)
    l_cache = nothing
  else
    ode_cache, l_cache = cache
  end

  ode_cache = update_cache!(ode_cache,op,tθ)

  afop = ThetaMethodAffineOperator(op,tθ,dtθ,u0,ode_cache)

  newmatrix = true
  l_cache = solve!(uf,solver.ls,afop,l_cache,newmatrix)

  if 0.0 < solver.θ < 1.0
    uf = uf*(1.0/solver.θ)-u0*((1-solver.θ)/solver.θ)
  end

  cache = (ode_cache, l_cache)

  tf = t0+dt
  return (uf,tf,cache)

end

function solve_step!(uf::AbstractVector,
                     solver::ThetaMethod,
                     op::ConstantODEOperator,
                     u0::AbstractVector,
                     t0::Real,
                     cache) # -> (uF,tF)

  dt = solver.dt
  solver.θ == 0.0 ? dtθ = dt : dtθ = dt*solver.θ
  tθ = t0+dtθ

  if cache === nothing
    ode_cache = allocate_cache(op)
    afop = ThetaMethodAffineOperator(op,tθ,dtθ,u0,ode_cache)
    l_cache = nothing
  else
    ode_cache, vθ, afop, l_cache = cache
  end

  ode_cache = update_cache!(ode_cache,op,tθ)

  residual!(afop.vector,odeop,tθ,u0,u0,ode_cache)

  newmatrix = false
  l_cache = solve!(uf,solver.ls,afop,l_cache,newmatrix)

  if 0.0 < solver.θ < 1.0
    uf = uf*(1.0/solver.θ)-u0*((1-solver.θ)/solver.θ)
  end

  cache = (ode_cache, afop, l_cache)

  tf = t0+dt
  return (uf,tf,cache)

end

"""
Affine operator that represents the θ-method affine operator at a
given time step, i.e., M(t)(u_n+θ-u_n)/dt + K(t)u_n+θ + b(t)
"""
function ThetaMethodAffineOperator(op::AffineODEOperator,tθ::Float64,dtθ::Float64,
                                   u0::AbstractVector,ode_cache)
  b = allocate_residual(odeop,u0)
  A = allocate_jacobian(odeop,u0)
  residual!(b,odeop,tθ,u0,v0,ode_cache)
  z = zero(eltype(A))
  fill_entries!(A,z)
  jacobian!(A,odeop,tθ,u0,u0,ode_cache)
  jacobian_t!(A,odeop,tθ,u0,u0,(1/op.dtθ),ode_cache)
  # santiagobadia : Not sure about the sign
  afop = AffineOperator(A,b)
end

function update_rhs!(b,op,tθ,ode_cache)
  residual!(b,odeop,tθ,u0,u0,ode_cache)
end
