function time_derivative(f::Function)
  _f(x) = t -> f(x,t)
  ∂tu(t) = x -> ForwardDiff.derivative(_f(x),t)
  ∂tu(x,t) = ∂tu(t)(x)
  ∂tu
end

const ∂t = time_derivative
