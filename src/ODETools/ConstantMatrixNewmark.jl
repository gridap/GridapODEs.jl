function solve_step!(
  x1::NTuple{3,AbstractVector},
  solver::Newmark,
  op::ConstantMatrixODEOperator,
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
    newmatrix = true

    # Allocate caches
    newmark_cache = allocate_cache(op,v0,a0)
    (v,a, ode_cache) = newmark_cache

    # Define Newmark operator
    newmark_affOp = NewmarkAffineOperator(op,t1,dt,γ,β,(u0,v0,a0),newmark_cache)

    # Allocate matrices and vectors
    A, b = _allocate_matrix_and_vector(op,x0,ode_cache)
    jacobian!(A,newmark_affOp,u1)

    # Create affine operator cache
    affOp_cache = (A,b,nothing)
  else
    newmatrix = true
    newmark_cache, affOp_cache = cache
  end

  # Unpack and update caches
  (v,a, ode_cache) = newmark_cache
  ode_cache = update_cache!(ode_cache,op,t1)
  A,b,l_cache = affOp_cache

  # Fill vector
  _vector!(b,newmark_affOp,u1)

  # Create affine operator with updated RHS
  affOp = AffineOperator(A,b)
  l_cache = solve!(u1,solver.nls,affOp,l_cache,newmatrix)

  # Update auxiliar variables
  u1 = u1 + u0
  v1 = γ/(β*dt)*(u1-u0) + (1-γ/β)*v0 + dt*(1-γ/(2*β))*a0
  a1 = 1.0/(β*dt^2)*(u1-u0) - 1.0/(β*dt)*v0 - (1-2*β)/(2*β)*a0

  # Pack caches
  affOp_cache = A,b,l_cache
  cache = (newmark_cache, affOp_cache)
  x1 = (u1,v1,a1)

  return (x1,t1,cache)

end
