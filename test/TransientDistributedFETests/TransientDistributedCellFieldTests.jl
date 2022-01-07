module TransientDistributedCellFieldTests

using Gridap
using PartitionedArrays
using GridapODEs.TransientDistributedFETools
using GridapDistributed: DistributedCellField
using Test

parts = get_part_ids(sequential,(2,))
domain = (0,1)
cells = (4,)
𝒯 = CartesianDiscreteModel(parts,domain,cells)
Ω = Interior(𝒯)
dΩ = Measure(Ω,2)

f(t) = t^2
df(t) = 2t
ddf(t) = 2

a(t) = CellField(f(t),Ω)
da(t) = CellField(df(t),Ω)
dda(t) = CellField(ddf(t),Ω)
@test isa(a(0),DistributedCellField)
@test isa(da(0),DistributedCellField)
@test isa(dda(0),DistributedCellField)

b(t) = TransientCellField(a(t),(da(t),dda(t)))
@test isa(b(0),TransientDistributedCellField)
@test isa(b(0),TransientSingleFieldDistributedCellField)

db(t) = ∂t(b(t))
@test isa(db(0),TransientDistributedCellField)
@test isa(db(0),TransientSingleFieldDistributedCellField)

ddb(t) = ∂t(db(t))
@test isa(ddb(0),TransientDistributedCellField)
@test isa(ddb(0),TransientSingleFieldDistributedCellField)

@test (∑(∫(a(0.5))dΩ)) ≈ 0.25
@test (∑(∫(da(0.5))dΩ)) ≈ 1.0
@test (∑(∫(dda(0.5))dΩ)) ≈ 2.0
@test (∑(∫(b(0.5))dΩ)) ≈ 0.25
@test (∑(∫(db(0.5))dΩ)) ≈ 1.0
@test (∑(∫(ddb(0.5))dΩ)) ≈ 2.0
@test (∑(∫(∂tt(b(0.5)))dΩ)) ≈ 2.0

end
