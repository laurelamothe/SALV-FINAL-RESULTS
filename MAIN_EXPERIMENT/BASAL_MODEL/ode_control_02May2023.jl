## Laure SALINITY MAIN EXP.
# Ostreo growth
# David/Laure - 02 May 2023

function ODE_ostreo(u,p,t)
    # Variables and parameters to be estimated
    r, K, ω  = p
    # ODE
    return r*u*(1-u/K) - ω*u
  end
  