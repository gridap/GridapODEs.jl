function solve_step!(
  x1::NTuple{3,AbstractVector},
  solver::Newmark,
  op::ConstantODEOperator,
  x0::NTuple{3,AbstractVector},
  t0::Real,
  cache) # -> (uF,tF)

  dt = solver.dt
  γ = solver.γ
  β = solver.β
  t1 = t0+dt
  u0, v0, a0 = x0
  u1, v1, a1 = x1

  if cache === nothing
    newmark_cache = allocate_cache(op,v0,a0)
    A, b = _allocate_matrix_and_vector(op,x0)
    b1 = similar(b)
    b1 .= 0.0
    M = allocate_jacobian(op,x0)
    newmatrix = true
    newmark_affOp = NewmarkAffineOperator(op,t1,dt,γ,β,(u0,v0,a0),newmark_cache, affOp_cache)
    _matrix_and_vector!(newmark_affOp,u1)
    A,b,l_cache = newmark_affOp.affop_cache
    M = _mass_matrix!(M,newmark_affOp,u1)
    affOp_cache = (A,b,M,)
  else
    newmark_cache, affOp_cache = cache
    newmatrix = false
  end

  (v,a, ode_cache) = newmark_cache
  ode_cache = update_cache!(ode_cache,op,t1)
  A,b,M,l_cache = affOp_cache
  b1 = b + M*u0
  affOp = AffineOperator(A,b1)
  l_cache = solve!(u1,solver.nls,affOp,l_cache,newmatrix)

  v1 = γ/(β*dt)*(u1-u0) + (1-γ/β)*v0 + dt*(1-γ/(2*β))*a0
  a1 = 1.0/(β*dt^2)*(u1-u0) - 1.0/(β*dt)*v0 - (1-2*β)/(2*β)*a0

  affOp_cache = A,b,M,l_cache
  cache = (newmark_cache, affOp_cache)
  x1 = (u1,v1,a1)

  return (x1,t1,cache)

end
