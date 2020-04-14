#@fverdugo this is the part that needs some improvements.
# - Most of the comments in other files are related with the design of this part.
#   Specially, with the cache objects.
# - We also want Dirichlet conditions represented by several objects.
#   For the moment, it is assumed that the Dirichlet conditions are described by a single function.
# - Some extra considerations have to be taken into account for MultiField problems.
#
# We can define a transient trial fe space as a FE space whose underlying dirichlet values 
# depend on time, but all other behavior is stationary (e.g., shape functions, cell dof basis, etc.)
# This is enough for the method of lines.
#
# Specifically, this is the behavior we need for a transient trial fe space.
#
# - evaluation at a given time t (We have this now,
#   but perhaps we need to modify a bit the API to avoid side effects
#   when evaluating the same object in two different time points, see below.)
#
# - Differentiation with respect to t (We have this now)
#
# - In may context we need a FE Space but we don't care about
# the underlying Dirichlet data  (e.g., for building an assembler, for getting
#  the shape functions when evaluating the jacobian in the TransientFEOperator, etc.).
#  This is the missing part that would help to improve other parts of the code.
#  We can implement (U::TransientTrialFESpace)(::Nothing) to be called as U(nothing)
#  that returns a space with arbitrary Dirichlet values to be used in contexts when we do not
#  care about the specific Dirichlet data. We can also implement
#  HomogeneousTrialFESpace(::TransientTrialFESpace) which could also be used in similar contexts.
#

struct TransientTrialFESpace
  space::FESpace
  dirichlet_t#::Vector{<:Function} #@fverdugo extend to dirichlet conditions represented by several objects. Not all of them need to be functions.
  dirichlet_values::AbstractVector #@fverdugo this is very dangerous since we cannot evaluate this
  # space at two different times simultaneously, which is a strong limitation.
  # We need to store this vector somewhere else (e.g., in the ode cache).
end

# function TransientTrialFESpace(U::FESpace,dir_fun)
  # TransientTrialFESpace(U,[dir_fun])
# end

function TransientTrialFESpace(U::FESpace,objects)
 dirichlet_vals = similar(get_dirichlet_values(U))
  TransientTrialFESpace(U,objects,dirichlet_vals)
end

function (fes::TransientTrialFESpace)(t::Real)
  # dir_t = [ud(t) for ud in fes.dirichlet_t]
  dir_t = fes.dirichlet_t(t)
  #@fverdugo This constructor needs to be called TrialFESpace! since the first argument is mutated.
  # To be changed in Gridap
  TrialFESpace(fes.dirichlet_values,fes.space,dir_t)
end

(tfes::FESpace)(t::Real) = tfes

∂t(fes::TransientTrialFESpace) = TransientTrialFESpace(fes.space,∂t.(fes.dirichlet_t))

∂t(fes::TrialFESpace) = HomogeneousTrialFESpace(fes.space)

#@fverdugo this has to be moved to Gridap
function HomogeneousTrialFESpace(U::FESpace)
  # @santiagobadia : To be improved
  TrialFESpace(U,0.0)
  # @fverdugo this should to be implemented as:
  #
  # dirichlet_values = zero_dirichlet_values(U)
  # TrialFESpace(dirichlet_values,U)
  #
  # The constructor TrialFESpace(dirichlet_values,U) is missing. To be added in Gridap.
  # In fact, this should be the inner constructor of the TrialFESpace struct.
end

function test_transient_trial_fe_space(Uh::TransientTrialFESpace)
  Uh0=Uh(0.0)
  Uh00=Uh0(0.0)
  Uht=∂t(Uh)
  Uht0=Uht(0.0)
  Uh0t = ∂t(Uh0)
  true
end
