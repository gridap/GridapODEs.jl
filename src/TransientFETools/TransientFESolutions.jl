# """
# A general TransientFESolution
# """
struct TransientFESolution
  solver::TransientFESolver
  op::TransientFEOperator
  u0#::CellField # What should I put here, SingleFieldFEFunction???
  # We need a name for all FEFunction, why not called FEFunction ???
  t0::Real
  tF::Real
end
