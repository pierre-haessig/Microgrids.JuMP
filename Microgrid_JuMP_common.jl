# Common code for Microgrid definition, simulation and optimization
# Pierre Haessig, September-October 2024
# This is a script and function-wrapped version of the code developped in Microgrid_optimization_JuMP.ipynb
# which is then reused in other notebooks

using Microgrids
using JuMP

println("Microgrid optimization with JuMP common functions:")
println("- CRF")
println("- ts_reduction")
println("- g_tan")
println("- cons_Xann_usage!")
println("- build_optim_mg_stage!")
println("- setup_optim_mg_jump")
println("- diagnostics_mg_jump")
println("- optim_mg_jump")
println("- simulate_alg")

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
CRF(0.05, 25)

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

"""
    build_optim_mg_stage!(mg::Microgrid, model_data::Dict{String,Any})

Build one stage of microgrid sizing optimization JuMP model

All named variables and constraints are unregistered from the JuMP model,
thus the use of the `model_data` Dict.

# Parameters:
- `mg`: base Microgrid description (with 1kW(h) ratings) for e.g. price parameters and load data
- `model_data`::Dict to store variable references and constraints.
  - a JuMP `Model` should be included in `model_data["model"]`
  - it should also contain the following model parameters (see `setup_optim_mg_jump`):
    `shed_max`, `ndays`, `fixed_lifetimes`, `gen_hours_assum`, `relax_gain` and `z_tan`

⚠ Sizing maximum bounds are taken as global variables:
- `power_rated_gen_max`
- `energy_rated_sto_max`
- `power_rated_pv_max`
- `power_rated_wind_max`

Also, a global `ts_reduction(x, ndays)` function is needed.
"""
function build_optim_mg_stage!(mg::Microgrid, model_data::Dict{String,Any})
    # Retrieve model parameters from the Dict:
    md = model_data # shortcut Dict name

    shed_max = md["shed_max"]
    ndays = md["ndays"]
    fixed_lifetimes = md["fixed_lifetimes"]
    gen_hours_assum = md["gen_hours_assum"]
    relax_gain = md["relax_gain"]
    z_tan = md["z_tan"]

    println("Building stage problem with $ndays days...")
    dt = mg.project.timestep
    discount_rate = mg.project.discount_rate
    CRFproj(T) = CRF(mg.project.discount_rate, T)

    K = ndays*24 # h
    ts_reduction_ndays(x) = ts_reduction(x, ndays)

    Pload = mg.load |> ts_reduction_ndays
    Eload_desired = sum(Pload)*dt*365/ndays
    md["Pload"] = Pload
    md["Eload_desired"] = Eload_desired

    # (works because the rated power in mg are set to 1 kW)
    cf_pv   = production(mg.nondispatchables[1]) |> ts_reduction_ndays
    cf_wind = production(mg.nondispatchables[2]) |> ts_reduction_ndays

    ### JuMP model definition
    model = md["model"] # JuMP model


    ##  Sizing variables
    md["power_rated_gen"]  = @variable(model, 0 <= power_rated_gen  <= power_rated_gen_max)
    md["energy_rated_sto"] = @variable(model, 0 <= energy_rated_sto <= energy_rated_sto_max)
    md["power_rated_pv"]   = @variable(model, 0 <= power_rated_pv   <= power_rated_pv_max)
    md["power_rated_wind"] = @variable(model, 0 <= power_rated_wind <= power_rated_gen_max)

    ## Power flows for each component

    # Non dispatchable sources (renewables)
    md["pv_potential"] = @variable(model, pv_potential[1:K])
    @constraint(model, pv_potential .== power_rated_pv*cf_pv)
    md["wind_potential"] = @variable(model, wind_potential[1:K])
    @constraint(model, wind_potential .== power_rated_wind*cf_wind);
    md["renew_potential"] = renew_potential = pv_potential + wind_potential

    # Desired net load (before spillage of excess renewables and load shedding)
    md["Pnl"] = @variable(model, Pnl[1:K])
    @constraint(model, Pnl .== Pload .- renew_potential)

    # Renewables spillage and load shedding
    md["Pspill"] = @variable(model, Pspill[1:K] >= 0)
    md["Pshed"]  = @variable(model, Pshed[1:K] >= 0)

    # Dispatchable generator
    md["Pgen"] = @variable(model, Pgen[1:K] >= 0)
    @constraint(model, Pgen .<= power_rated_gen)

    ## Energy storage
    # Charge and discharge power (net of losses)
    md["Psto_cha"] = @variable(model, Psto_cha[1:K] >= 0)
    md["Psto_dis"] = @variable(model, Psto_dis[1:K] >= 0)
    # Energy level and its bounds
    md["Esto"] = @variable(model, Esto[1:K+1])
    @constraint(model, Esto .<= energy_rated_sto) # SoCmax = 1 implicitly
    @constraint(model, Esto .>= mg.storage.SoC_min*energy_rated_sto) # often 0
    # Power constraints
    @constraint(model, Psto_cha .<= mg.storage.charge_rate    * energy_rated_sto) #
    @constraint(model, Psto_dis .<= mg.storage.discharge_rate * energy_rated_sto) # double the solving time with HiGHS and MOI 1.30.0, when there is the Pgen penalty!!
    # variant that embeds the two preceding ones, using that fact that Psto_cha(k) * Psto_dis(k) = 0
    @constraint(model, Psto_cha/mg.storage.charge_rate + Psto_dis/mg.storage.discharge_rate  .<=  energy_rated_sto)
    # Evolution of the State of Energy,
    # with piecewise linear in P loss term (aka constant efficiency model)
    a = mg.storage.loss_factor
    md["stodyn"] = @constraint(model, [k=1:K],
        Esto[k+1] == Esto[k] + (Psto_cha[k] - Psto_dis[k] - a*(Psto_cha[k]+Psto_dis[k]))*dt,
        base_name = "stodyn"
    )
    # Storage cyclicity
    @constraint(model, Esto[K+1] == Esto[1])
    # not implemented: force initial SoC # TODO: activate to make it comparable with BB
    #@constraint(model, Esto[1] == mg.storage.SoC_ini * energy_rated_sto)

    ## Power balance
    md["balance"] = @constraint(model, balance,
        Pgen + (Psto_dis - Psto_cha) - Pspill .== Pnl - Pshed,
    )

    if shed_max == 0.0
        println("zero load shedding allowed")
        fix.(Pshed, 0.0; force=true);
    else
        println("load shedding allowed: $shed_max")
        Eshed = sum(Pshed)*dt * 365/ndays
        @constraint(model, Eshed <= Eload_desired*shed_max)
    end

    ### Costs
    # Generator costs
    md["Pgen_rated_ann"] = @variable(model, Pgen_rated_ann >= 0) # annualized size
    md["Egen"] = Egen = sum(Pgen)*dt * 365/ndays # Generator yearly energy

    if fixed_lifetimes
        print("Fixed generator lifetime hypothesis: ")
        gen_lifetime = mg.generator.lifetime_hours / gen_hours_assum # years
        println("$gen_lifetime y, assuming $gen_hours_assum  h/y of usage")
        md["cons_Pgen_rated_ann_CRFgen_lifetime"] = @constraint(model,
            Pgen_rated_ann == power_rated_gen * CRFproj(gen_lifetime))
    else
        println("Usage-dependent generator lifetime model (relax_gain=", relax_gain,")")
        md["Ugen"] = @variable(model, Ugen >= 0) # cumulated usage
        println("Xann constraints with z_tan=$z_tan")
        @constraint(model, Ugen == Egen*relax_gain/mg.generator.lifetime_hours); # kW/y
        cpwl_gen = cons_Xann_usage!(model,
            Pgen_rated_ann, power_rated_gen, Ugen,
            discount_rate, z_tan)
    end
    # question for O&M wit fixed lifetime: use assumed fixed gen hours ?
    # or use the convexified hours (Egen*relax_gain) like for Ugen?
    md["Cgen_om"] = Cgen_om = fixed_lifetimes ?
        mg.generator.om_price_hours * gen_hours_assum * power_rated_gen :
        mg.generator.om_price_hours * relax_gain * Egen

    md["Cgen_fuel"] = Cgen_fuel = mg.generator.fuel_price * mg.generator.fuel_slope * Egen;# $/y
    md["Cgen"] = Cgen = mg.generator.investment_price * Pgen_rated_ann +
                        Cgen_om + Cgen_fuel # $/y

    # Battery costs
    md["Esto_rated_ann"] = @variable(model, Esto_rated_ann >= 0) # annualized size
    md["E_through_sto"] = E_through_sto = (sum(Psto_cha) + sum(Psto_dis))*dt * 365/ndays # cumulated throughput
    md["Csto"] = Csto = mg.storage.investment_price * Esto_rated_ann +
                        mg.storage.om_price * energy_rated_sto
    # A) Effect of calendar lifetime:
    CRFsto_cal = CRFproj(mg.storage.lifetime_calendar)
    md["cons_Esto_rated_ann_CRFsto_cal"] = @constraint(model,
        Esto_rated_ann >= energy_rated_sto*CRFsto_cal)
    # B) Effect of cycling
    if ~fixed_lifetimes
        md["Usto"] = @variable(model, Usto >= 0) # cumulated usage
        @constraint(model, Usto == E_through_sto/(2*mg.storage.lifetime_cycles))
        cpwl_sto = cons_Xann_usage!(model,
            Esto_rated_ann, energy_rated_sto, Usto,
            discount_rate, z_tan)
    end

    # Wind and solar costs
    pv = mg.nondispatchables[1]
    md["Cpv"] = Cpv = pv.investment_price * power_rated_pv * CRFproj(pv.lifetime) +
                      pv.om_price * power_rated_pv
    wind = mg.nondispatchables[2]
    md["Cwind"] = Cwind = wind.investment_price * power_rated_wind * CRFproj(wind.lifetime) +
                          wind.om_price * power_rated_wind

    # Total cost
    md["Cann"] = Cann = Cgen + Csto +  Cpv + Cwind
    md["LCOE"] = Cann/Eload_desired

    # Unregister all variables
    for var in keys(object_dictionary(model))
        unregister(model, var)
    end

    return nothing
end

"""Setup microgrid sizing optimization model using JuMP

# Parameters:
- `optimizer`: JuMP compatible optimization solver

Keyword parameter:
- `create_mg_base` function returning a `Microgrid` project description
  (e.g. for load, components price and lifetime)

Optional keyword parameters:
- `shed_max`: load shedding upper bound, as a fraction of desired load energy, in [0,1] (default to 0.0)
- `ndays`: time series reduction (default to 365: no data reducation)
- `fixed_lifetimes`: true or false (default)
- `gen_hours_assum`: assumed generator operation hours when `fixed_lifetimes` is `true`,
   in 0 – 8760 h/y (default 2000 h/y)
- `relax_gain`: 1.0 by default. Increase to try to compensation the underestimation
  of e.g. generator operating hours due to linearization.
- `z_tan`: tangent points for Xann PWL approx, see `cons_Xann_usage!` function

(all these parameters gets stored in the `model_data` output Dict)

⚠ a global `create_microgrid(x; create_mg_base)` function is needed. TO BE REMOVED.

Create a JuMP model, populates it with `build_optim_mg_stage!` and returns:
- mg: base microgrid
- model_data Dict, with optimization variables and named constraints
"""
function setup_optim_mg_jump(optimizer; create_mg_base,
        shed_max=0.0,
        ndays=365,
        fixed_lifetimes=false,
        gen_hours_assum = 2000.,
        relax_gain = 1.0,
        z_tan = [0.20, 0.28, 0.37, 0.50, 0.68, 1.0, 1.7, 4.0],
    )

    # base Microgrid with 1kW(h) ratings
    mg = create_microgrid([1., 1., 1., 1.]; create_mg_base)

    model_data = Dict{String,Any}()

    # JuMP model setup
    model = Model(optimizer)
    set_silent(model)
    model_data["model"] = model

    # Store model parameters in the Dict:
    model_data["shed_max"] = shed_max
    model_data["ndays"] = ndays
    model_data["fixed_lifetimes"] = fixed_lifetimes
    model_data["gen_hours_assum"] = gen_hours_assum
    model_data["relax_gain"] = relax_gain
    model_data["z_tan"] = z_tan

    # Populate the model
    build_optim_mg_stage!(mg, model_data)

    # Set optimization objective
    @objective(model, Min, model_data["Cann"])

    return mg, model_data
end

"""Diagnostics about microgrid optimization with JuMP

About operation stats and economic performance, for generator and storage.
e.g. lifetime and CRF (including the effect of annualized sized pwl approximation),

Returns a hierarchical NamedTuple
"""
function diagnostics_mg_jump(mg, model_data, ndays, relax_gain)
    md = model_data
    dt = mg.project.timestep
    Eload_desired = md["Eload_desired"]
    CRFproj(T) = CRF(mg.project.discount_rate, T)
    # Generator diagnostics
    Pgen = value.(md["Pgen"])
    Egen = value(md["Egen"])
    power_rated_gen = value(md["power_rated_gen"])
    # operation hours/y
    gen_hours = sum(Pgen .> 1e-6*power_rated_gen)*dt*365/ndays
    gen_hours_lin = relax_gain*Egen/power_rated_gen
    gen_lifetime = mg.generator.lifetime_hours / gen_hours
    gen_lifetime_hlin = mg.generator.lifetime_hours / gen_hours_lin
    if "Ugen" in keys(md)
        @assert isapprox(gen_lifetime_hlin, power_rated_gen/value(md["Ugen"]); rtol=1e-8)
    end
    # Storage diagnostics
    energy_rated_sto = value(md["energy_rated_sto"])
    E_through_sto = value(md["E_through_sto"])
    sto_cycles = value(E_through_sto/(2*energy_rated_sto)) # c/year
    sto_lifetime_cycles = mg.storage.lifetime_cycles / sto_cycles
    if "Usto" in keys(md)
        @assert isapprox(sto_lifetime_cycles, energy_rated_sto/value(md["Usto"]); rtol=1e-8)
    end
    sto_lifetime = min(sto_lifetime_cycles, mg.storage.lifetime_calendar)

    diagnostics = (
        generator = (
            cost_share = value(md["Cgen"]/md["Cann"]),
            cost_share_fuel = value(md["Cgen_fuel"]/md["Cann"]),
            energy = Egen,
            load_share = Egen/Eload_desired,
            hours = gen_hours,
            hours_lin = gen_hours_lin,
            lifetime = gen_lifetime,
            lifetime_hlin = gen_lifetime_hlin,
            CRF = CRFproj(gen_lifetime),
            CRF_hlin = CRFproj(gen_lifetime_hlin),
            CRF_hlin_pwl = value(md["Pgen_rated_ann"] / power_rated_gen)
        ),
        storage = (
            cost_share = value(md["Csto"]/md["Cann"]),
            energy_through = E_through_sto,
            load_share = 0.5*E_through_sto/Eload_desired,
            cycles = sto_cycles,
            lifetime_cycles = sto_lifetime_cycles,
            lifetime = sto_lifetime,
            CRF = CRFproj(sto_lifetime),
            CRF_pwl = value(md["Esto_rated_ann"]/energy_rated_sto)
        )
    )
    return diagnostics
end

"""Optimize sizing of microgrid using JuMP

and extract results

# Parameters
- `optimizer`: JuMP compatible optimization solver

Keyword parameter:
- `create_mg_base` function returning a `Microgrid` project description
  (e.g. for load, components price and lifetime)

Optional keyword parameter:
- optional keyword parameters of the `setup_optim_mg_jump` function
- `model_custom`: a function taking `model_data` and `mg_base` as arguments,
  which can modify the model before its optimization

Returns:
xopt, LCOE_opt, diagnostics, traj, model_data
"""
function optim_mg_jump(optimizer; create_mg_base,
        shed_max=0.0,
        ndays=365,
        fixed_lifetimes=false,
        gen_hours_assum = 2000.,
        relax_gain = 1.0,
        z_tan = [0.20, 0.28, 0.37, 0.50, 0.68, 1.0, 1.7, 4.0],
        model_custom=nothing
    )

    # Setup optimization model
    mg, model_data = setup_optim_mg_jump(optimizer;
        shed_max, ndays, fixed_lifetimes, gen_hours_assum, relax_gain, z_tan,
        create_mg_base
    )
    md = model_data
    model = md["model"]
    LCOE = md["LCOE"]

    # Optional model customization
    if model_custom !== nothing
        model_custom(model_data, mg)
    end

    # Run optimization
    @time JuMP.optimize!(model)

    # Extract results
    LCOE_opt = value(LCOE)
    xopt = value.([
        md["power_rated_gen"]
        md["energy_rated_sto"]
        md["power_rated_pv"]
        md["power_rated_wind"]
    ])

    # Paste back sizing in Microgrid project description
    mg.generator.power_rated = value(md["power_rated_gen"])
    mg.storage.energy_rated = value(md["energy_rated_sto"])
    mg.nondispatchables[1].power_rated = value(md["power_rated_pv"])
    mg.nondispatchables[2].power_rated = value(md["power_rated_wind"])

    diagnostics = diagnostics_mg_jump(mg, model_data, ndays, relax_gain)

    traj = (
        Pgen = value.(md["Pgen"]),
        Psto = value.(md["Psto_dis"] - md["Psto_cha"]),
        Esto = value.(md["Esto"]),
        Pshed = value.(md["Pshed"])
    )
    return xopt, LCOE_opt, diagnostics, traj, md, mg
end

"""
    simulate_alg(mg::Microgrid, md, smoothing::Smoothing=NoSmoothing)

simulate solution of Algebraic model (stored as JuMP variables in `model_data`)
using Microgrids's  simulator.

To get comparable results, `smoothing=Smoothing(transition=1.0, gain=relax_gain)`
should be used.

This replicates `Microgrids.simulate`, except that the `operation` step is bypassed.
"""
function simulate_alg(mg::Microgrid, md, smoothing::Smoothing=NoSmoothing)
    # Collect microgrid operation trajectories
    Pnl = value.(md["Pnl"])
    Pshed = value.(md["Pshed"])
    renew_potential = value.(md["renew_potential"])
    Pgen = value.(md["Pgen"])
    Esto = value.(md["Esto"])
    Psto = value.(md["Psto_dis"] - md["Psto_cha"])
    Pspill = value.(md["Pspill"])
    Psto_dmax = zero(Psto)
    Psto_cmax = zero(Psto)

    oper_traj = OperationTraj(Pnl, Pshed, renew_potential, Pgen, Esto, Psto, Psto_dmax, Psto_cmax, Pspill)

    # Aggregate the operation variables
    oper_stats = aggregation(mg, oper_traj, smoothing)

    # Eval the microgrid costs
    mg_costs = economics(mg, oper_stats)

    return (traj=oper_traj, stats=oper_stats, costs=mg_costs)
end
