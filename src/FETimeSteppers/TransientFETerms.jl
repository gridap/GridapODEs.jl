struct TransientFETerm <: FETerm # assuming first order for now
  res::Function
  jac_unk::Function
  jac_unk_t::Function
  trian::Triangulation
  quad::CellQuadrature
end
# @santiagobadia : Not all the terms must be transient in the TransientFEOperator, but it would
# imply to create jacobian_unk_t methods that return 0 for these other terms
# It can be done in GridapTimeStepper as an extension of the FEOperator
# interface in Gridap

get_cell_jacobian_unk = get_cell_jacobian

get_cell_jacobian_unk_t(::FETerm,uh,uh_t,du_t,v) = @notimplemented #or 0
