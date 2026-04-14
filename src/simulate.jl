
using Microgrids


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
    LoH   = zeros(length(Pnl))
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