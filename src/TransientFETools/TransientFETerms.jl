# @santiagobadia : Create an AbstractTransientFETerm (?)

struct TransientFETerm <: FETerm # assuming first order for now
  res::Function
  jac::Function
  jac_t::Function
  trian::Triangulation
  quad::CellQuadrature
end
# @santiagobadia : If we use subtyping, it could lead to errors since a standard
# FEOperatorFromTerms could take a transient term, nonsense ???

# @santiagobadia : Not all the terms must be transient in the TransientFEOperator, but it would
# imply to create jacobian_unk_t methods that return 0 for these other terms
# It can be done in GridapTimeStepper as an extension of the FEOperator
# interface in Gridap

get_cell_jacobian_t(::FETerm,uh,uh_t,du_t,v) = @notimplemented #or 0
# @santiagobadia :  Should this term depend on t, idem for _jacobian
