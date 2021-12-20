struct TransientMultiFieldCellField{A} <: TransientCellField
  transient_single_fields::Vector{<:TransientCellField} # used to iterate
  cellfield::A
  derivatives::Tuple
end

MultiFieldTypes = Union{MultiFieldCellField,MultiFieldFEFunction}

function TransientCellField(multi_field::MultiFieldTypes,derivatives::Tuple)
  transient_single_fields = TransientCellField[]
  for ifield in 1:num_fields(multi_field)
    single_field = multi_field[ifield]
    single_derivatives = ()
    for ifield_derivatives in derivatives
      single_derivatives = (single_derivatives...,getindex(ifield_derivatives,ifield))
    end
    transient_single_field = TransientSingleFieldCellField(single_field,single_derivatives)
    push!(transient_single_fields,transient_single_field)
  end
  TransientMultiFieldCellField(transient_single_fields,multi_field,derivatives)
end

function get_data(f::TransientMultiFieldCellField)
  s = """
  Function get_data is not implemented for TransientMultiFieldCellField at this moment.
  You need to extract the individual fields and then evaluate them separatelly.

  If ever implement this, evaluating a `MultiFieldCellField` directly would provide,
  at each evaluation point, a tuple with the value of the different fields.
  """
  @notimplemented s
end

function get_triangulation(f::TransientMultiFieldCellField)
  s1 = first(f.single_fields)
  trian = get_triangulation(s1)
  @check all(map(i->trian===get_triangulation(i),f.single_fields))
  trian
end
DomainStyle(::Type{TransientMultiFieldCellField{DS}}) where DS = DS()
num_fields(a::TransientMultiFieldCellField) = length(a.transient_single_fields)
Base.getindex(a::TransientMultiFieldCellField,i::Integer) = a.transient_single_fields[i]
Base.iterate(a::TransientMultiFieldCellField)  = iterate(a.transient_single_fields)
Base.iterate(a::TransientMultiFieldCellField,state)  = iterate(a.transient_single_fields,state)

# Time derivative
function ∂t(f::TransientMultiFieldCellField)
  transient_single_field_derivatives = TransientCellField[]
  for transient_single_field in f.transient_single_fields
    push!(transient_single_field_derivatives,∂t(transient_single_field))
  end
  cellfield, derivatives = first_and_tail(f.derivatives)
  TransientMultiFieldCellField(transient_single_field_derivatives,cellfield,derivatives)
end

∂tt(f::TransientMultiFieldCellField) = ∂t(∂t(f))

function Base.view(a::TransientMultiFieldCellField,indices::Vector{<:Int})
  transient_single_fields = TransientCellField[]
  cellfield = MultiFieldCellField(a.cellfield[indices],DomainStyle(a.cellfield))
  derivatives = ()
  for derivative in a.derivatives
    derivatives = (derivatives...,MultiFieldCellField(derivative[indices],DomainStyle(derivative)))
  end
  for ifield in indices
    single_field = a.cellfield[ifield]
    single_derivatives = ()
    for ifield_derivatives in a.derivatives
      single_derivatives = (single_derivatives...,getindex(ifield_derivatives,ifield))
    end
    transient_single_field = TransientSingleFieldCellField(single_field,single_derivatives)
    push!(transient_single_fields,transient_single_field)
  end
  TransientMultiFieldCellField(transient_single_fields,cellfield,derivatives)
end
