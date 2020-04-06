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

function FETerm(
  res::Function, jac::Function, jac_t::Function, trian::Triangulation, quad::CellQuadrature)
  TransientFETermFromIntegration(res,jac,jac_t,trian,quad)
end

function get_cell_residual(tr::TransientFETermFromIntegration,t::Real,uh,uh_t,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,tr.trian)
  _uh = restrict(uh,tr.trian)
  _uh_t = restrict(uh_t,tr.trian)
  integrate(tr.res(_uh,_uh_t,_v),t.trian,t.quad)
end

function get_cell_jacobian(tr::TransientFETermFromIntegration,t::Real,uh,uh_t,du,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  _v = restrict(v,tr.trian)
  _uh = restrict(uh,tr.trian)
  _uh_t = restrict(uh_t,tr.trian)
  _du = restrict(du,tr.trian)
  integrate(tr.jac(t,_uh,_uh_t,_du,_v),t.trian,t.quad)
end

function get_cell_jacobian_t(tr::TransientFETermFromIntegration,t::Real,uh,uh_t,du_t,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du_t)
  _v = restrict(v,tr.trian)
  _uh = restrict(uh,tr.trian)
  _uh_t = restrict(uh_t,tr.trian)
  _du_t = restrict(du_t,tr.trian)
  integrate(tr.jac_t(t,_uh,_uh_t,_du_t,_v),tr.trian,tr.quad)
end

function get_cell_values(tr::TransientFETermFromIntegration,uhd)
  @assert is_a_fe_function(uhd)
  cellvals = get_cell_values(uhd)
  reindex(cellvals,tr.trian)
end

function get_cell_id(t::TransientFETermFromIntegration)
  get_cell_id(t.trian)
end
