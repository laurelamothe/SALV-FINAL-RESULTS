constant(c, Sal) = c .* (Sal .* 0.0 .+ 1)

linear(a, b, Sal) = a .* Sal .+ b

poly(a, b, c, Sal) = a .* Sal .^ 2 .+ b .* Sal .+ c

logistic(a, k, s0, Sal) = k ./ (1 .+ exp.(.-a .* (Sal .- s0)))

normal(a, σ, µ, Sal) = a .* (1 ./ (σ .* sqrt(2 * pi)) .* exp.(-1 / 2 .* ((Sal .- µ) ./ σ) .^ 2))


hyperfunctions = Dict(
    "r" => poly,
    "Kr" => linear,
    "ϕ" => linear,
    "η" => normal,
    "K" => normal,
    "µ" => constant,
    "ϵ" => linear,
    "β" => constant,
    "ρ" => constant,
    "Δ" => constant)

