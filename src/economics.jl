using Microgrids



"""Capital recovery factor for discount rate `i` and duration `T`
CRF is such that Cann = NPC*CRF
"""
function CRF(i,T)
    if i != 0.0
        a = (1+i)^T
        return i*a/(a-1)
    else
        return 1/T
    end
end

"""intercept and slope of tangent to function g:z → 1/(1-e^(-1/z)) at point `z0>0`

Usage: with gi, g1 = g_tan(z0),
then g(z) ~ gi + g1*z around z0
gi: intercept at z=0
g1: slope
"""
function g_tan(z0)
    z=z0
    emiz = exp(-1/z)
    g = 1/(1-emiz) # g(z)
    g1 = emiz/(z^2*(1-emiz)^2) # g'(z)
    ginter = g - g1*z
    return ginter, g1
end

"""Piecewise linear approximation of “annualized size” function of size and usage

Parameters:
- Xann, X, U: JuMP variables or expressions
- discount_rate: Real number in [0,1] (negative discount not implemented)
- z_tan: list of z0 values where to use tangents of g(z)

Recommended values:
- z_tan = [0.5]: up to 5% CRF underestimation error
- z_tan = [0.28, 0.5, 1.0]: 1.6% error
- z_tan = [0.20, 0.28, 0.37, 0.50, 0.68, 1.0, 1.7, 4.0]: 0.23% error

Tangents at z=0 and z→+∞ are also added, so no need to add these.

Returns the vector of constraints
"""
function cons_Xann_usage!(model, Xann, X, U, discount_rate, z_tan=[0.28, 0.5, 1.0])
    r = discount_rate/log(1+discount_rate)
    # stacked tangents: (ginter, g1) for each z0 in z_tan, including z=0 and +infty
    tan_params = [
        (1.0, 0.) # z0 = 0, yields Xann >= X*discount_rate
        [(g_tan(z0)[1], g_tan(z0)[2]) for z0 in z_tan]
        (0.5, 1.0) # z0 -> +inftyn yields Xann >= X*discount_rate*0.5 + U*r
    ]
    ntan = length(tan_params)
    cvec = @constraint(model, [k=1:ntan],
        Xann >= X*discount_rate*tan_params[k][1]
              + U*r*tan_params[k][2]
    )
    cname = "c_" * name(Xann) *"_"* name(X) *"_"* name(U)
    set_name.(cvec, [cname * "[$k]" for k=1:ntan])
    return cvec
end


"""Calculation of the real cost of this model after optimization

Parameters:
- mg : Microgrid optimized
- md : Model with power rated and storage
- td : time reduced to calculate overall quantity of H2

Returns the size of the storage needed, the cost that comes with it, and the amount of H2 over time
"""


function Q_hydro_overtime(mg::Microgrid, md, td::Vector{Float64} = zeros(365*24))
    #Variables needed
    Pele = value.(md["Pele"])
    cons_rate_elyz = mg.electrolyzer[1].consumption_slope
    Pfc   = value.(md["Pfc"])
    cons_rate_fc = mg.dispatchables.fuel_cell[1].consumption_slope
    dt = td[2] - td[1]
    K = length(td)
    hytank= mg.tanks.h2Tank

    Q = cumsum(((Pele[1:K] ./ cons_rate_elyz) .- (Pfc[1:K] .* cons_rate_fc)) .* dt)
    range_Q = maximum(Q) - minimum(Q)
    Q_om_price = hytank.capacity * hytank.om_price

    cost_Q = range_Q * hytank.investment_price * CRF(mg.project.discount_rate,hytank.lifetime) + Q_om_price
    return (Q, range_Q, cost_Q)
end