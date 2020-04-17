
"""
A single field FE space with transient Dirichlet data (see Multifield below).
"""
struct TransientTrialFESpace
  space::SingleFieldFESpace
  dirichlet_t::Union{Function,Vector{<:Function}}
  Ud0::TrialFESpace

  function TransientTrialFESpace(space::SingleFieldFESpace,dirichlet_t::Union{Function,Vector{<:Function}})
    Ud0 = HomogeneousTrialFESpace(space)
    new(space,dirichlet_t,Ud0)
  end
end

"""
Time evaluation without allocating Dirichlet vals
"""
function evaluate!(Ut::TrialFESpace,U::TransientTrialFESpace,t::Real)
  if isa(U.dirichlet_t,Vector)
    objects_at_t = map( o->o(t), U.dirichlet_t)
  else
    objects_at_t = U.dirichlet_t(t)
  end
  TrialFESpace!(Ut,objects_at_t)
  Ut
end

"""
Allocate the space to be used as first argument in evaluate!
"""
function allocate_trial_space(U::TransientTrialFESpace)
    HomogeneousTrialFESpace(U.space)
end

"""
Time evaluation allocating Dirichlet vals
"""
function evaluate(U::TransientTrialFESpace,t::Real)
  Ut = allocate_trial_space(U)
  evaluate!(Ut,U,t)
  Ut
end

"""
We can evaluate at `nothing` when we do not care about the Dirichlet vals
"""
function evaluate(U::TransientTrialFESpace,t::Nothing)
  U.Ud0
end

"""
Functor-like evaluation. It allocates Dirichlet vals in general.
"""
(U::TransientTrialFESpace)(t) = evaluate(U,t)

"""
Time derivative of the Dirichlet functions
"""
∂t(U::TransientTrialFESpace) = TransientTrialFESpace(U.space,∂t.(U.dirichlet_t))


# Testing the interface

function test_transient_trial_fe_space(Uh)
  UhX = Uh(nothing)
  @test isa(UhX,FESpace)
  Uh0 = allocate_trial_space(Uh)
  Uh0 = evaluate!(Uh0,Uh,0.0)
  @test isa(Uh0,FESpace)
  Uh0 = evaluate(Uh,0.0)
  @test isa(Uh0,FESpace)
  Uh0 = Uh(0.0)
  @test isa(Uh0,FESpace)
  Uht=∂t(Uh)
  Uht0=Uht(0.0)
  @test isa(Uht0,FESpace)
  true
end

# Define the TransientTrialFESpace interface for stationary spaces

function evaluate!(Ut::FESpace,U::FESpace,t::Real)
  U
end

function allocate_trial_space(U::FESpace)
  U
end

function evaluate(U::FESpace,t::Real)
  U
end

function evaluate(U::FESpace,t::Nothing)
  U
end

(U::FESpace)(t) = U

∂t(U::FESpace) = HomogeneousTrialFESpace(U)

# Define the interface for MultiField

struct TransientMultiFieldTrialFESpace
  spaces::Vector
end

function MultiFieldFESpace(spaces::Vector)
  TransientMultiFieldTrialFESpace(spaces)
end

function evaluate!(Ut::MultiFieldFESpace,U::TransientMultiFieldTrialFESpace,t::Real)
  spaces_at_t = [evaluate!(Ut.spaces[i],U.spaces[i],t) for i in 1:length(U.spaces)]
  MultiFieldFESpace(spaces_at_t)
end

function allocate_trial_space(U::TransientMultiFieldTrialFESpace)
  spaces = allocate_trial_space.(U.spaces)
  MultiFieldFESpace(spaces)
end

function evaluate(U::TransientMultiFieldTrialFESpace,t::Real)
  Ut = allocate_trial_space(U)
  evaluate!(Ut,U,t)
  Ut
end

(U::TransientMultiFieldTrialFESpace)(t) = evaluate(U,t)

function ∂t(U::TransientMultiFieldTrialFESpace)
  spaces = ∂t.(U.spaces)
  MultiFieldFESpace(spaces)
end
