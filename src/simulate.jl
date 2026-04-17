
using Microgrids


"""reduce the year time series `x` to `ndays` ≤ 365
sr=24: data sampling rate / day. 24 means hourly time step.

The method is *basic*:
it samples `ndays` days evenly spaced by (365 ÷ `ndays`) days

with ndays=365, returns the original series
"""
function ts_reduction(x, ndays, sr=24)
    out = zeros(ndays*sr)
    @assert ndays<=365
    Δdays = 365 ÷ ndays
    for i in 1:ndays
        offset_in = (i-1)*Δdays*sr
        offset_out = (i-1)*sr
        out[offset_out+1:offset_out+sr] = x[offset_in+1:offset_in+sr]
    end
    return out
end

"""
    simulate_alg(mg::Microgrid, md, smoothing::Smoothing=NoSmoothing)

simulate solution of Algebraic model (stored as JuMP variables in `model_data`)
using Microgrids's  simulator.

Discontinuous statistics can optionally be relaxed (smoothed)
using the relaxation parameter `ε`:
- 0.0 means no relaxation (default value)
- 1.0 yields the strongest relaxation

Using relaxation (`ε` > 0) is recommended when using gradient-based optimization
and then a “small enough” value between 0.05 and 0.30 is suggested.

This replicates `Microgrids.simulate`, except that the `operation` step is bypassed.
"""
function simulate_alg(mg::Microgrid, md, ε::Real=0.0)

    Pnl   = value.(md["Pnl"])
    Pshed = value.(md["Pshed"])
    Prenew = value.(md["renew_potential"])
    Pgen = zeros(length(Pnl))
    Pfc   = value.(md["Pfc"])
    Pele = value.(md["Pele"])
    Phb   = zeros(length(Pnl))
    Esto = value.(md["Esto"])
    Psto = value.(md["Psto_dis"] .- md["Psto_cha"])
    LoH   = value.(md["LoH"])
    LoF   = zeros(length(Pnl))
    Pcurt = value.(md["Pspill"])
    Pdump = zeros(length(Pnl))

    oper_traj = OperationTraj(Pnl, Pshed, Prenew, Pgen, Pfc, Pele,
    Phb, Esto, Psto, LoH, LoF, Pcurt, Pdump
    )

    oper_stats = aggregation(mg, oper_traj, ε)
    mg_costs = economics(mg, oper_stats)

    return (traj=oper_traj, stats=oper_stats, costs=mg_costs)
end