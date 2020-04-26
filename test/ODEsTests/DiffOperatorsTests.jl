module DiffOperatorsTests

using Gridap
using Test
using GridapODEs.ODETools: ∂t

using ForwardDiff

f(x,t) = 5*x[1]*x[2]+x[2]^2*t^3
∂tf(x,t) = x[2]^2*3*t^2

tv = rand(Float64)
xv = Point(rand(Float64,2)...)
@test ∂tf(xv,tv) ≈ ∂t(f)(xv,tv)


F(x,t) = [5*x[1]*x[2],x[2]^2*t^3]
∂tF(x,t) = [0.0,x[2]^2*3*t^2]

tv = rand(Float64)
xv = Point(rand(Float64,2)...)
@test ∂tF(xv,tv) ≈ ∂t(F)(xv,tv)

end #module
# @santiagobadia : Trying to compute time derivatives for MultiValue
# without success below
# ##
# X = Gridap.Inference
# z = X.return_type(Ft(tv),typeof(xv))
#
# # ForwardDiff for multivariate functions
#
# # The case of univariate function is clear, as already implemented in Gridap
# h(x) = VectorValue(5*x[1]*x[2],x[2]^2)
#
# dh(x) = ForwardDiff.jacobian(y->h(y).array,x.array)
# dh(x) = ForwardDiff.jacobian(y->h(y).array,x.array)
# dh(xv)
#
# # As soon as H does not depend on t it is OK
# H(x,t) = VectorValue(5*x[1]*x[2],x[2]^2)
# Ht(t) = x -> H(x,t)
# Hx(x) = t -> H(x,t)
# dH(t) = x -> ForwardDiff.derivative(y->Hx(x)(y).array,t)
# dH(tv)(xv)
#
# # But what I would expect to work for multivariate functions does not
# H(x,t) = VectorValue(5*x[1]*x[2],x[2]^2+t)
# Ht(t) = x -> H(x,t)
# Hx(x) = t -> H(x,t)
# dH(t) = x -> ForwardDiff.derivative(y->Hx(x)(y).array,t)
# dH(t) = x -> ForwardDiff.derivative(Hx(x),t)
# dH(tv)(xv)
#
# # Whereas everything smooth without VectorValue, i.e., arrays
# K(x,t) =    [5*x[1]*x[2],x[2]^2+t]
# Kx(x) = t -> K(x,t)
# dK(t) = x -> ForwardDiff.derivative(y->Kx(x)(y),t)
# dK(tv)(xv)
#
# # The VectorValue treatment is OK
# A(t) = VectorValue(5*t,2*t)
# ForwardDiff.derivative(y->A(y).array,tv)
#
# # As above, no problem for practically univariate functions
# B(x,t) = A(t)
# ForwardDiff.derivative(y->B(xv,y).array,tv)
#
# # It also works (sure)
# C(x) = t -> A(t)
# ForwardDiff.derivative(y->C(xv)(y).array,tv)
#
# # But it does not
# D(x,t) = VectorValue(5*x[1]*x[2],x[2]^2+t)
# ForwardDiff.derivative(y->D(xv,y).array,tv)
#
# # Neither this (as expected)
# E(x) = t -> VectorValue(5*x[1]*x[2],x[2]^2+t)
# ForwardDiff.derivative(y->E(xv)(y).array,tv)
#
# # That is surprising
# xv = Point(rand(Float64,2)...)
# FF = t -> VectorValue(5*xv[1]*xv[2],xv[2]^2+t)
# ForwardDiff.derivative(y->FF(y).array,tv)
#
# # Or even in a function (sth about the global scope...)
# function mytry(xv,D)
#     FF(t) = D(xv,t)
#     ForwardDiff.derivative(y->FF(y).array,tv)
# end
# GG = mytry(xv,D)
#
# ##
#
# FV(x,t) = VectorValue(5*x[1]*x[2],x[2]^2*t^3)
# ∂tFV(x,t) = VectorValue(0.0,x[2]^2*3*t^2)
#
# FVx(x) = t -> FV(x,t)
# FVt(t) = x -> FV(x,t)
#
# dtFVx(t) = x -> ForwardDiff.derivative(y->FVx(x)(y).array,t)
# # dtFVx(t) = VectorValue(x -> ForwardDiff.derivative(FVx(x),t))
# dtFV(x,t) = dtFVx(t)(x)
#
# tv = rand(Float64)
# xv = Point(rand(Float64,2)...)
# # @test ∂tFV(xv,tv) ≈
# dtFVx(tv)(xv)
# dH(tv)(xv)
#
#
# ##
# _time_derivative_f(f,x,t,z) = ForwardDiff.derivative(f(x),t)
# _time_derivative_f(f,x,t,z::VectorValue) = VectorValue(ForwardDiff.derivative(f(x),t))
#
# X = Gridap.Inference
# time_derivative(f)
#
# time_derivative_f(f,x,t)
#
# z = zero(X.return_type(ft(tv),typeof(xv)))
# @test _time_derivative_f(fx,xv,tv,z) ≈ ∂tf(xv,tv)
#
# z = zero(X.return_type(FVt(tv),typeof(xv)))
# @test
#
# VectorValue( x -> ForwardDiff.derivative(FVx(xv),tv))
#
#
# _time_derivative_f(FVx,xv,tv,z)
# ∂tFV(xv,tv)
#
