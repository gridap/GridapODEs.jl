# struct TransientMassAffineFETermFromIntegration <: TransientFETermFromIntegration
#   m::Function
#   a::Function
#   b::Function
#   trian::Triangulation
#   quad::CellQuadrature
# end
#
# function TransientAffineFETerm(
#   m::Function, a::Function, b::Function, trian::Triangulation, quad::CellQuadrature)
#   TransientMassAffineFETermFromIntegration(m,a,b,trian,quad)
# end
#
# function get_cell_residual(tr::TransientMassAffineFETermFromIntegration,t::Real,uh,uh_t,v)
#   @assert is_a_fe_function(uh)
#   @assert is_a_fe_function(uh_t)
#   @assert is_a_fe_cell_basis(v)
#   _v = restrict(v,tr.trian)
#   _uh = restrict(uh,tr.trian)
#   _uh_t = restrict(uh_t,tr.trian)
#   integrate(tr.m(t,_uh_t,_v)+tr.a(t,_uh,_v)-tr.b(t,_v),tr.trian,tr.quad)
# end
#
# function get_cell_jacobian(tr::TransientMassAffineFETermFromIntegration,t::Real,uh,uh_t,du,v)
#   @assert is_a_fe_function(uh)
#   @assert is_a_fe_function(uh_t)
#   @assert is_a_fe_cell_basis(v)
#   @assert is_a_fe_cell_basis(du)
#   _v = restrict(v,tr.trian)
#   _uh = restrict(uh,tr.trian)
#   _uh_t = restrict(uh_t,tr.trian)
#   _du = restrict(du,tr.trian)
#   integrate(tr.a(t,_du,_v),tr.trian,tr.quad)
# end
#
# function get_cell_jacobian_t(tr::TransientMassAffineFETermFromIntegration,t::Real,uh,uh_t,du_t,v,duht_du::Real)
#   @assert is_a_fe_function(uh)
#   @assert is_a_fe_function(uh_t)
#   @assert is_a_fe_cell_basis(v)
#   @assert is_a_fe_cell_basis(du_t)
#   _v = restrict(v,tr.trian)
#   _uh = restrict(uh,tr.trian)
#   _uh_t = restrict(uh_t,tr.trian)
#   _du_t = restrict(du_t,tr.trian)
#   integrate(duht_du*tr.m(t,_du_t,_v),tr.trian,tr.quad)
# end
#
# struct TransientFEMassFromIntegration <: TransientMassAffineFETermFromIntegration
#   m::Function
#   trian::Triangulation
#   quad::CellQuadrature
# end
#
# # santiagobadia : Without Transient no dispatching possible
# function TransienFEMassTerm(
#   m::Function, trian::Triangulation, quad::CellQuadrature)
#   TransientFEMassFromIntegration(m,trian,quad)
# end
#
# function get_cell_residual(tr::TransientFEMassFromIntegration,t::Real,uh,uh_t,v)
#   @assert is_a_fe_function(uh)
#   @assert is_a_fe_function(uh_t)
#   @assert is_a_fe_cell_basis(v)
#   _v = restrict(v,tr.trian)
#   _uh = restrict(uh,tr.trian)
#   _uh_t = restrict(uh_t,tr.trian)
#   integrate(tr.m(t,_uh_t,_v),tr.trian,tr.quad)
# end
#
# function get_cell_jacobian(tr::TransientFEMassFromIntegration,t::Real,uh,uh_t,du,v)
#   nothing
# end
#
# struct TransientAffineFETermFromIntegration <: TransientMassAffineFETermFromIntegration
#   a::Function
#   b::Function
#   trian::Triangulation
#   quad::CellQuadrature
# end
#
# # santiagobadia : Without Transient no dispatching possible
# function TransientAffineFETerm(
#   a::Function, b::Function, trian::Triangulation, quad::CellQuadrature)
#   TransientAffineFETermFromIntegration(a,b,trian,quad)
# end
#
# function get_cell_residual(tr::TransientAffineFETermFromIntegration,t::Real,uh,uh_t,v)
#   @assert is_a_fe_function(uh)
#   @assert is_a_fe_function(uh_t)
#   @assert is_a_fe_cell_basis(v)
#   _v = restrict(v,tr.trian)
#   _uh = restrict(uh,tr.trian)
#   _uh_t = restrict(uh_t,tr.trian)
#   integrate(tr.a(t,_uh,_v)-tr.b(t,_v),tr.trian,tr.quad)
# end
#
# function get_cell_jacobian(tr::TransientAffineFETermFromIntegration,t::Real,uh,uh_t,du,v)
#   @assert is_a_fe_function(uh)
#   @assert is_a_fe_function(uh_t)
#   @assert is_a_fe_cell_basis(v)
#   @assert is_a_fe_cell_basis(du)
#   _v = restrict(v,tr.trian)
#   _uh = restrict(uh,tr.trian)
#   _uh_t = restrict(uh_t,tr.trian)
#   _du = restrict(du,tr.trian)
#   integrate(tr.a(t,_du,_v),tr.trian,tr.quad)
# end
#
# function get_cell_jacobian_t(tr::TransientAffineFETermFromIntegration,t::Real,uh,uh_t,du_t,v,duht_du::Real)
#   nothing
# end
#
# struct TransientLinearFETermFromIntegration <: TransientAffineFETermFromIntegration
#   a::Function
#   trian::Triangulation
#   quad::CellQuadrature
# end
#
# # santiagobadia : Without Transient no dispatching possible
# function TransientLinearFETerm(
#   a::Function, trian::Triangulation, quad::CellQuadrature)
#   TransientLinearFETermFromIntegration(a,b,trian,quad)
# end
#
# function get_cell_residual(tr::TransientLinearFETermFromIntegration,t::Real,uh,uh_t,v)
#   @assert is_a_fe_function(uh)
#   @assert is_a_fe_function(uh_t)
#   @assert is_a_fe_cell_basis(v)
#   _v = restrict(v,tr.trian)
#   _uh = restrict(uh,tr.trian)
#   _uh_t = restrict(uh_t,tr.trian)
#   integrate(tr.a(t,_uh,_v),tr.trian,tr.quad)
# end
#
# struct TransientFESourceFromIntegration <: TransientAffineFETermFromIntegration
#   b::Function
#   trian::Triangulation
#   quad::CellQuadrature
# end
#
# # santiagobadia : Without Transient no dispatching possible
# function TransientFESource(
#   b::Function, trian::Triangulation, quad::CellQuadrature)
#   TransientFESourceFromIntegration(b,trian,quad)
# end
#
# function get_cell_residual(tr::TransientFESourceFromIntegration,t::Real,uh,uh_t,v)
#   @assert is_a_fe_function(uh)
#   @assert is_a_fe_function(uh_t)
#   @assert is_a_fe_cell_basis(v)
#   _v = restrict(v,tr.trian)
#   integrate(-tr.b(t,_v),tr.trian,tr.quad)
# end
#
# function get_cell_jacobian(tr::TransientFESourceFromIntegration,t::Real,uh,uh_t,du,v)
#   nothing
# end
