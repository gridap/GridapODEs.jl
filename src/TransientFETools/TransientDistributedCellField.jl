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

# Differential Operations
gradient(f::TransientDistributedCellField) = gradient(f.cellfield)
∇∇(f::TransientDistributedCellField) = ∇∇(f.cellfield)

# Unary ops
for op in (:symmetric_part,:inv,:det,:abs,:abs2,:+,:-,:tr,:transpose,:adjoint,:grad2curl,:real,:imag,:conj)
  @eval begin
    ($op)(a::TransientDistributedCellField) = Operation($op)(a.cellfield)
  end
end

# Binary ops
for op in (:inner,:outer,:double_contraction,:+,:-,:*,:cross,:dot,:/)
  @eval begin
    ($op)(a::TransientDistributedCellField,b::TransientDistributedCellField) = Operation($op)(a.cellfield,b.cellfield)
    ($op)(a::TransientDistributedCellField,b::DistributedCellField) = Operation($op)(a.cellfield,b)
    ($op)(a::DistributedCellField,b::TransientDistributedCellField) = Operation($op)(a,b.cellfield)
    ($op)(a::TransientDistributedCellField,b::Number) = Operation($op)(a.cellfield,b)
    ($op)(a::Number,b::TransientDistributedCellField) = Operation($op)(a,b.cellfield)
    ($op)(a::TransientDistributedCellField,b::Function) = Operation($op)(a.cellfield,b)
    ($op)(a::Function,b::TransientDistributedCellField) = Operation($op)(a,b.cellfield)
  end
end

Base.broadcasted(f,a::TransientDistributedCellField,b::TransientDistributedCellField) = Operation((i,j)->f.(i,j))(a.cellfield,b.cellfield)
Base.broadcasted(f,a::TransientDistributedCellField,b::DistributedCellField) = Operation((i,j)->f.(i,j))(a.cellfield,b)
Base.broadcasted(f,a::DistributedCellField,b::TransientDistributedCellField) = Operation((i,j)->f.(i,j))(a,b.cellfield)
Base.broadcasted(f,a::Number,b::TransientDistributedCellField) = Operation((i,j)->f.(i,j))(a,b.cellfield)
Base.broadcasted(f,a::TransientDistributedCellField,b::Number) = Operation((i,j)->f.(i,j))(a.cellfield,b)
Base.broadcasted(f,a::Function,b::TransientDistributedCellField) = Operation((i,j)->f.(i,j))(a,b.cellfield)
Base.broadcasted(f,a::TransientDistributedCellField,b::Function) = Operation((i,j)->f.(i,j))(a.cellfield,b)
Base.broadcasted(::typeof(*),::typeof(∇),f::TransientDistributedCellField) = Operation(Fields._extract_grad_diag)(∇(f.cellfield))
Base.broadcasted(::typeof(*),s::Fields.ShiftedNabla,f::TransientDistributedCellField) = Operation(Fields._extract_grad_diag)(s(f.cellfield))
