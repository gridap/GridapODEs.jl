abstract type TransientFETerm <: FETerm end

get_cell_jacobian_t(::FETerm,uh,uh_t,du_t,v) = @notimplemented #or 0

get_cell_jacobian_t(::TransientFETerm,uh,uh_t,du_t,v) = @notimplemented #or 0

struct TransientFETermFromIntegration <: TransientFETerm
  res::Function
  jac::Function
  jac_t::Function
  trian::Triangulation
  quad::CellQuadrature
end

get_cell_jacobian_t(::TransientFETermFromIntegration,uh,uh_t,du_t,v) = @notimplemented
