# A concrete version of ODEOperator from a FEOperator
struct ODEOperatorFromFEOperator <: ODEOperator
  op::FEOperator
end
