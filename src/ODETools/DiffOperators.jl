function time_derivative(f::Function)
  fx(x) = t -> f(x,t)
  ∂tf(t) = x -> ForwardDiff.derivative(fx(x),t)
  ∂tf(x,t) = ∂tf(t)(x)
  ∂tf
end

const ∂t = time_derivative

# @santiagobadia: Trying to create time derivatives for
# MultiValue without success
function array_time_derivative(f::Function)
  function tdt_f(f,xv,tv)
    ty = return_type(f(tv),typeof(xv))
    _time_derivative(f,xv,tv,ty)
  end
end

function _time_derivative(f::Function,x,t,ty::Float64)
  fx(t) = f(xv,t)
  ∂tf = ForwardDiff.derivative(fx,t)
  # ∂tf(x,t) = ∂tf(t)(x)
  ∂tf
end

function _time_derivative(f::Function,::Type{Vector})
  fx(x) = t -> f(x,t)
  ∂tf(t) = x -> VectorValue(ForwardDiff.derivative(fx(x),t)...)
  ∂tf(x,t) = ∂tf(t)(x)
end
