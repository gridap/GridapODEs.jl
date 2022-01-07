# Transient Distributed CellField
abstract type TransientDistributedCellField <: DistributedCellDatum  end

# Transient SingleField
struct TransientSingleFieldDistributedCellField{A} <: TransientDistributedCellField
  cellfield::A
  derivatives::Tuple
end

DistributedSingleFieldTypes = Union{DistributedCellField,DistributedSingleFieldFEFunction}

# Constructors
function TransientCellField(single_field::DistributedSingleFieldTypes,derivatives::Tuple)
  TransientSingleFieldDistributedCellField(single_field,derivatives)
end

# Time derivative
function ∂t(f::TransientDistributedCellField)
  cellfield, derivatives = first_and_tail(f.derivatives)
  TransientCellField(cellfield,derivatives)
end

∂tt(f::TransientDistributedCellField) = ∂t(∂t(f))

# Integration related
function Fields.integrate(f::TransientDistributedCellField,b::DistributedMeasure)
  integrate(f.cellfield,b)
end
