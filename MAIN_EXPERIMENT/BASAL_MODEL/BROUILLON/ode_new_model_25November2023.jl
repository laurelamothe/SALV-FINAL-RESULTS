## Laure SALINITY MAIN EXP.
# Ostreo + OtVs system
# David/Laure - 04 May 2023

function ODE_system(du, u, p, t)
    # Variables and parameters to be estimated
    S, R, I, V = u
    N = S + R + I
    r, μ, K, e, α, ϵ, η, β, δ = p
    # ODE
    du[1] = (1 - μ) * r * S * (1 - N / K) - e * α * r * S * V # Susceptible
    du[2] = r * (μ * S + ϵ * R) * (1 - N / K)          # Resistant
    du[3] = e * α * r * S * V - η * r * I   # Infected
    du[4] = β * η * r * I - δ * V - e * S * V # VLPs
end


