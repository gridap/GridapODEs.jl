# Transient CellField
struct TransientCellField{A} <: CellField
  cellfield::A
  derivatives::Tuple{Vararg{A}}
end

# CellField methods
get_data(f::TransientCellField) = get_data(f.cellfield)
get_triangulation(f::TransientCellField) = get_triangulation(f.cellfield)
DomainStyle(::Type{<:TransientCellField{A}}) where A = DomainStyle(A)
gradient(f::TransientCellField) = gradient(f.cellfield)
∇∇(f::TransientCellField) = ∇∇(f.cellfield)
change_domain(f::TransientCellField,trian::Triangulation,target_domain::DomainStyle) = change_domain(f.cellfield,trian,target_domain)

# MultiFieldCellField methods
num_fields(f::TransientCellField{A}) where A = num_fields(f.cellfield)
function Base.getindex(f::TransientCellField{A},i::Integer) where A 
  singleCellfield = getindex(f.single_cellfields,i)
  singleDerivatives = (getindex(i_derivatives,i) for i_derivatives in f.derivatives)
  TransientCellField(singleCellfield,singleDerivatives)
end
function Base.iterate(f::TransientCellField{A}) where A
  cellField, state1 = iterate(f.cellfield)
  derivatives = ()
  state2 = []
  for i_derivatives in f.derivatives
    single_derivative, single_state = iterate(i_derivatives)
    derivatives = (derivatives...,single_derivative)
    push!(state2,single_state)
  end
  TransientCellField(cellField,derivatives),(state1,state2)
end
function Base.iterate(f::TransientCellField{A},state) where A
  state1, state2 = state
  cellField, state1 = iterate(f.cellfield,state[1])
  derivatives = ()
  for (i,i_derivatives) in enumerate(f.derivatives)
    single_derivative, state2[i] = iterate(i_derivatives)
    derivatives = (derivatives...,single_derivative)
  end
  TransientCellField(cellField,derivatives),(state1, state2)
end

# Transient FEBasis
struct TransientFEBasis{A} <: FEBasis
  febasis::A
  derivatives::Tuple{Vararg{A}}
end

# FEBasis methods
get_data(f::TransientFEBasis) = get_data(f.febasis)
get_triangulation(f::TransientFEBasis) = get_triangulation(f.febasis)
DomainStyle(::Type{<:TransientFEBasis{A}}) where A = DomainStyle(A)
BasisStyle(::Type{<:TransientFEBasis{A}}) where A = BasisStyle(A)
gradient(f::TransientFEBasis) = gradient(f.febasis)
∇∇(f::TransientFEBasis) = ∇∇(f.febasis)
change_domain(f::TransientFEBasis,trian::Triangulation,target_domain::DomainStyle) = change_domain(f.febasis,trian,target_domain)

# Time derivative
function ∂t(f::Union{TransientCellField,TransientFEBasis})
  cellfield, derivatives = first_and_tail(f.derivatives)
  TransientCellField(cellfield,derivatives)
end

∂tt(f::Union{TransientCellField,TransientFEBasis}) = ∂t(∂t(f::Union{TransientCellField,TransientFEBasis}))