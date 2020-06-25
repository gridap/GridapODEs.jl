struct ConstantMassAffineFETermFromIntegration <: TransientFETermFromIntegration
  m::Function
  a::Function
  b::Function
  trian::Triangulation
  quad::CellQuadrature
end

function MassAffineFETerm(
  m::Function, a::Function, b::Function, trian::Triangulation, quad::CellQuadrature)
  MassAffineFETermFromIntegration(m,a,b,trian,quad)
end

function get_cell_residual(tr::TransientMassAffineFETermFromIntegration,t::Real,uh,uh_t,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,tr.trian)
  _uh = restrict(uh,tr.trian)
  _uh_t = restrict(uh_t,tr.trian)
  integrate(tr.m(_uh_t,_v)+tr.a(_uh,_v)-tr.b(_v),tr.trian,tr.quad)
end

function get_cell_jacobian(tr::TransientMassAffineFETermFromIntegration,t::Real,uh,uh_t,du,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du)
  _v = restrict(v,tr.trian)
  _uh = restrict(uh,tr.trian)
  _uh_t = restrict(uh_t,tr.trian)
  _du = restrict(du,tr.trian)
  integrate(tr.a(_du,_v),tr.trian,tr.quad)
end

function get_cell_jacobian_t(tr::TransientMassAffineFETermFromIntegration,t::Real,uh,uh_t,du_t,v,duht_du::Real)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(du_t)
  _v = restrict(v,tr.trian)
  _uh = restrict(uh,tr.trian)
  _uh_t = restrict(uh_t,tr.trian)
  _du_t = restrict(du_t,tr.trian)
  integrate(duht_du*tr.m(_du_t,_v),tr.trian,tr.quad)
end

function get_cell_matrix_and_vector(t::AffineFETermFromIntegration,uhd,u,v)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  cellmat = get_cell_matrix(t,u,v)
  cellvec = _get_cell_vector_tmp_hack(cellmat,t,v,uhd) #TODO
  cellvals = get_cell_values(t,uhd)
  _setup_cell_matrix_and_vector(cellmat,cellvec,cellvals)
end

function get_cell_matrix(t::AffineFETermFromIntegration,u,v)
  @assert is_a_fe_cell_basis(v)
  @assert is_a_fe_cell_basis(u)
  _v = restrict(v,t.trian)
  _u = restrict(u,t.trian)
  integrate(t.biform(_u,_v),t.trian,t.quad)
end

function get_cell_vector(t::AffineFETermFromIntegration,uhd,v)
  @assert is_a_fe_function(uhd)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  _uhd = restrict(uhd,t.trian)
  integrate(t.liform(_v)-t.biform(_uhd,_v),t.trian,t.quad)
end

function get_cell_vector(t::AffineFETermFromIntegration,v)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,t.trian)
  integrate(t.liform(_v),t.trian,t.quad)
end

function get_cell_id(t::AffineFETermFromIntegration)
  get_cell_id(t.trian)
end

function get_cell_values(t::AffineFETermFromIntegration,uhd)
  @assert is_a_fe_function(uhd)
  cellvals = get_cell_values(uhd)
  reindex(cellvals,t.trian)
end

struct ConstantFEMassFromIntegration <: ConstantMassAffineFETermFromIntegration
  m::Function
  trian::Triangulation
  quad::CellQuadrature
end

# santiagobadia : Without Transient no dispatching possible
function TransienFEMassTerm(
  m::Function, trian::Triangulation, quad::CellQuadrature)
  TransientFEMassFromIntegration(m,trian,quad)
end

function get_cell_residual(tr::TransientFEMassFromIntegration,t::Real,uh,uh_t,v)
  @assert is_a_fe_function(uh)
  @assert is_a_fe_function(uh_t)
  @assert is_a_fe_cell_basis(v)
  _v = restrict(v,tr.trian)
  _uh = restrict(uh,tr.trian)
  _uh_t = restrict(uh_t,tr.trian)
  integrate(tr.m(t,_uh_t,_v),tr.trian,tr.quad)
end

function get_cell_jacobian(tr::TransientFEMassFromIntegration,t::Real,uh,uh_t,du,v)
  nothing
end
