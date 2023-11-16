## Laure SALINITY MAIN EXP.
# Ostreo + OtVs system
# David/Laure - 04 May 2023

function ODE_system(du, u, p, t)
    # Variables and parameters to be estimated
    S, R, I, V = u
    N = S + R + I
    r, μ, Ks, Kr, ω, ϵ, ϕ, η, β, δ = p
    # ODE
    du[1] = (1 - μ) * r * S * (1 - N / Ks) - ω * S - ϕ * S * V # Susceptible
    #du[2] = r * (μ * S + ϵ * R) * (1 - N / Kr) - ω * R         # Resistant
    du[2] = r * μ * S * (1 - N / Ks) + r * ϵ * R * (1 - N / Kr) - ω * R         # Resistant
    du[3] = ϕ * S * V - ω * I - η * I   # Infected
    du[4] = β * η * I - δ * V - ϕ * S * V # VLPs
end


