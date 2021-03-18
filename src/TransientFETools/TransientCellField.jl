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