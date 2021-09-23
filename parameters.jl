function parameters(
    ;
    dt::Rational = 1//24,
    T::Float64 = 30.,
    c::Float64 = .6,
    u1::Float64 = 2., # 10cm/s ~ 10km/d ~ 1-10 hab widths per day
    k1::Float64 = 4., # 1000m^2/s ~ 100 km^2/d ~ 1-100 hab widths per day
    r::Float64 = .25
    )

    p = Dict()

    # time
    p[:T] = T
    p[:Tc] = c*T
    p[:dt] = dt
    p[:tvec] = float.(0:dt:T)
    p[:NT] = length(p[:tvec])

    # mortality
    p[:μ_lo] = .0125 # per-day mortality rate
    p[:μ_hi] = .05

    # movement
    p[:u0] = -r*u1
    p[:u1] = u1
    p[:k0] = sqrt(r)*k1
    p[:k1] = k1
    p[:ℓ] = 1. # cbl width

    # passive movement
    p[:λ0] = 14.
    p[:λ1] = 1.

    return(p)
end
