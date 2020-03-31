# Do we need a TransientFEOperatorFromTerms or the only thing we need to do is
# to extend its interface. I would consider the second option

function jacobian_unk_t!(A::AbstractMatrix,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian_unk_t(uh,uh_t,du_t,v,op.terms)
  assemble_matrix!(A,op.assem, cellmats, cellidsrows, cellidscols)
  A
end
