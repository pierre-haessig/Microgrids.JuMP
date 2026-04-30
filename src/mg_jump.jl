using JuMP
using HiGHS
using Microgrids
include("../data/Microgrid_model_data.jl")
include("economics.jl")
include("simulate.jl")
include("create_microgrid.jl")



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
    `shed_max`, `ndays`, `fixed_lifetimes`, `fc_hours_assum`, `relax_gain` and `z_tan`

⚠ Sizing maximum bounds are taken as global variables:
- `power_rated_fc_max`
- `power_rated_ele`
- `energy_rated_sto_max`
- `power_rated_pv_max`
- `power_rated_wind_max`

Also, a global `ts_reduction(x, ndays)` function is needed.
"""

function build_optim_mg_stage!(mg::Microgrid, model_data::Dict{String,Any}, H2_price = 0.0, infinite_storage::Bool = false)
    md = model_data
    shed_max = md["shed_max"]
    ndays = md["ndays"]
    fixed_lifetimes = md["fixed_lifetimes"]
    fc_hours_assum = md["fc_hours_assum"]
    relax_gain = md["relax_gain"]
    z_tan = md["z_tan"]

    dt = mg.project.timestep
    discount_rate = mg.project.discount_rate
    CRFproj(T) = CRF(mg.project.discount_rate, T)

    K = ndays*24
    ts_reduction_ndays(x) = ts_reduction(x, ndays)

    Pload = mg.load |> ts_reduction_ndays
    Eload_desired = sum(Pload)*dt*365/ndays
    md["Pload"] = Pload
    md["Eload_desired"] = Eload_desired

    raw_pv = production(mg.nondispatchables[1])
    cf_pv = (raw_pv ./ maximum(raw_pv)) |> ts_reduction_ndays

    raw_wind = production(mg.nondispatchables[2])
    cf_wind = (raw_wind ./ maximum(raw_wind)) |> ts_reduction_ndays

    model = md["model"]

    md["power_rated_fc"]  = @variable(model, 0 <= power_rated_fc)
    md["power_rated_ele"] = @variable(model, 0 <= power_rated_ele)
    md["energy_rated_sto"] = @variable(model, 0 <= energy_rated_sto)

   if infinite_storage
        Ctank = 0.0
    else
        Ctank = @variable(model, 0 <= Ctank)
        md["Ctank"] = Ctank
    end

    md["power_rated_pv"]   = @variable(model, 0 <= power_rated_pv)
    md["power_rated_wind"] = @variable(model, 0 <= power_rated_wind)

    md["pv_potential"] = @variable(model, pv_potential[1:K])
    @constraint(model, pv_potential .== power_rated_pv*cf_pv)
    md["wind_potential"] = @variable(model, wind_potential[1:K])
    @constraint(model, wind_potential .== power_rated_wind*cf_wind)
    md["renew_potential"] = renew_potential = pv_potential + wind_potential

    md["Pnl"] = @variable(model, Pnl[1:K])
    @constraint(model, Pnl .== Pload .- renew_potential)

    md["Pspill"] = @variable(model, Pspill[1:K] >= 0)
    md["Pshed"]  = @variable(model, Pshed[1:K] >= 0)

    md["Pfc"] = @variable(model, Pfc[1:K] >= 0)
    @constraint(model, Pfc .<= power_rated_fc)

    md["Pele"] = @variable(model, Pele[1:K] >= 0)
    @constraint(model, Pele .<= power_rated_ele)

    md["Psto_cha"] = @variable(model, Psto_cha[1:K] >= 0)
    md["Psto_dis"] = @variable(model, Psto_dis[1:K] >= 0)
    md["Esto"] = @variable(model, Esto[1:K+1])
    @constraint(model, Esto .<= energy_rated_sto)
    @constraint(model, Esto .>= mg.storage.SoC_min*energy_rated_sto)
    @constraint(model, Psto_cha .<= mg.storage.charge_rate    * energy_rated_sto)
    @constraint(model, Psto_dis .<= mg.storage.discharge_rate * energy_rated_sto)
    @constraint(model, Psto_cha/mg.storage.charge_rate + Psto_dis/mg.storage.discharge_rate  .<=  energy_rated_sto)

    a = mg.storage.loss_factor
    md["stodyn"] = @constraint(model, [k=1:K],
        Esto[k+1] == Esto[k] + (Psto_cha[k] - Psto_dis[k] - a*(Psto_cha[k]+Psto_dis[k]))*dt
    )
    @constraint(model, Esto[K+1] == Esto[1])

    if ~infinite_storage 
        hytank = mg.tanks.h2Tank
        fc = mg.dispatchables.fuel_cell[1]
        elec = mg.electrolyzer[1]

        md["LoH"] = @variable(model, LoH[1:K+1] >= 0)

        @constraint(model, LoH .<= hytank.max_filling_ratio * Ctank)
        @constraint(model, LoH .>= hytank.min_filling_ratio * Ctank)

        @constraint(model, [k=1:K],
            LoH[k+1] == LoH[k]
                - Pfc[k] * fc.consumption_slope*dt
                + (Pele[k]/elec.consumption_slope)*dt
        )
        @constraint(model, LoH[K+1] == LoH[1])
    end

    md["balance"] = @constraint(model,
        Pfc - Pele + (Psto_dis - Psto_cha) - Pspill .== Pnl - Pshed
    )

    if shed_max == 0.0
        fix.(Pshed, 0.0; force=true)
    else
        Eshed = sum(Pshed)*dt * 365/ndays
        @constraint(model, Eshed <= Eload_desired*shed_max)
    end

    md["Pfc_rated_ann"] = @variable(model, Pfc_rated_ann >= 0)
    md["Efc"] = Efc = sum(Pfc)*dt * 365/ndays

    if fixed_lifetimes
        fc_lifetime = mg.dispatchables.fuel_cell[1].lifetime_hours / fc_hours_assum
        @constraint(model, Pfc_rated_ann == power_rated_fc * CRFproj(fc_lifetime))
    else
        md["Ufc"] = @variable(model, Ufc >= 0)
        @constraint(model, Ufc == Efc*relax_gain/mg.dispatchables.fuel_cell[1].lifetime_hours)
        cons_Xann_usage!(model, Pfc_rated_ann, power_rated_fc, Ufc, discount_rate, z_tan)
    end


    md["Cfc_om"] = fixed_lifetimes ?
        mg.dispatchables.fuel_cell[1].om_price_hourly * fc_hours_assum * power_rated_fc :
        mg.dispatchables.fuel_cell[1].om_price_hourly * relax_gain * Efc

    md["Cfc_combustible"] = infinite_storage ? H2_price * mg.dispatchables.fuel_cell[1].consumption_slope * Efc : 0.0
    md["Cfc"] = mg.dispatchables.fuel_cell[1].investment_price * Pfc_rated_ann +
                md["Cfc_om"] + md["Cfc_combustible"]

    md["Pele_rated_ann"] = @variable(model, Pele_rated_ann >= 0)
    md["Eele"] = Eele = sum(Pele)*dt * 365/ndays

    if fixed_lifetimes
        ele_lifetime = mg.electrolyzer[1].lifetime_hours / fc_hours_assum
        @constraint(model, Pele_rated_ann == power_rated_ele * CRFproj(ele_lifetime))
    else
        md["Uele"] = @variable(model, Uele >= 0)
        @constraint(model, Uele == Eele*relax_gain/mg.electrolyzer[1].lifetime_hours)
        cons_Xann_usage!(model, Pele_rated_ann, power_rated_ele, Uele, discount_rate, z_tan)
    end

    md["Cele_om"] = fixed_lifetimes ?
        mg.electrolyzer[1].om_price_hourly * fc_hours_assum * power_rated_ele :
        mg.electrolyzer[1].om_price_hourly * relax_gain * Eele

    md["Cele_combustible"] = infinite_storage ? H2_price / mg.electrolyzer[1].consumption_slope * Eele : 0.0
    md["Cele"] = infinite_storage ?
        mg.electrolyzer[1].investment_price * Pele_rated_ann + md["Cele_om"] - md["Cele_combustible"] :
        mg.electrolyzer[1].investment_price * Pele_rated_ann + md["Cele_om"]

    md["Esto_rated_ann"] = @variable(model, Esto_rated_ann >= 0)
    md["E_through_sto"] = E_through_sto = (sum(Psto_cha) + sum(Psto_dis))*dt * 365/ndays

    md["Csto"] = mg.storage.investment_price * Esto_rated_ann +
                 mg.storage.om_price * energy_rated_sto

    CRFsto_cal = CRFproj(mg.storage.lifetime_calendar)
    @constraint(model, Esto_rated_ann >= energy_rated_sto*CRFsto_cal)

    if ~fixed_lifetimes
        md["Usto"] = @variable(model, Usto >= 0)
        @constraint(model, Usto == E_through_sto/(2*mg.storage.lifetime_cycles))
        cons_Xann_usage!(model, Esto_rated_ann, energy_rated_sto, Usto, discount_rate, z_tan)
    end

    hytank = mg.tanks.h2Tank
    if infinite_storage
        md["CH2"] = 0.0
    else
        md["CH2"] = hytank.investment_price*md["Ctank"]*CRFproj(hytank.lifetime) +
                    hytank.om_price*md["Ctank"]
    end

    pv = mg.nondispatchables[1]
    md["Cpv"] = pv.investment_price * power_rated_pv * CRFproj(pv.lifetime) +
                pv.om_price * power_rated_pv

    wind = mg.nondispatchables[2]
    md["Cwind"] = wind.investment_price * power_rated_wind * CRFproj(wind.lifetime) +
                  wind.om_price * power_rated_wind

    md["Cann"] = infinite_storage ?
        md["Cfc"] + md["Csto"] + md["Cpv"] + md["Cwind"] + md["Cele"] :
        md["Cfc"] + md["Csto"] + md["Cpv"] + md["Cwind"] + md["Cele"] + md["CH2"]

    md["LCOE"] = md["Cann"]/Eload_desired

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
- `fc_hours_assum`: assumed fuel_cell operation hours when `fixed_lifetimes` is `true`,
   in 0 – 8760 h/y (default 2000 h/y)
- `relax_gain`: 1.0 by default. Increase to try to compensation the underestimation
  of e.g. fuel_cell operating hours due to linearization.
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
        fc_hours_assum = 2000.,
        relax_gain = 1.0,
        z_tan = [0.20, 0.28, 0.37, 0.50, 0.68, 1.0, 1.7, 4.0],
        H2_price,
        infinite_storage,
    )

    mg = create_microgrid([1., 1., 1., 1., 1., 1., 1., 1.]; create_mg_base)

    model_data = Dict{String,Any}()

    model = Model(optimizer)
    set_silent(model)
    model_data["model"] = model

    model_data["shed_max"] = shed_max
    model_data["ndays"] = ndays
    model_data["fixed_lifetimes"] = fixed_lifetimes
    model_data["fc_hours_assum"] = fc_hours_assum
    model_data["relax_gain"] = relax_gain
    model_data["z_tan"] = z_tan

    build_optim_mg_stage!(mg, model_data, H2_price, infinite_storage)

    @objective(model, Min, model_data["Cann"])

    return mg, model_data
end

"""Diagnostics about microgrid optimization with JuMP

About operation stats and economic performance, for fuel_cell and storage.
e.g. lifetime and CRF (including the effect of annualized sized pwl approximation),

Returns a hierarchical NamedTuple
"""
function diagnostics_mg_jump(mg, model_data, ndays, relax_gain)
    md = model_data
    dt = mg.project.timestep
    Eload_desired = md["Eload_desired"]
    CRFproj(T) = CRF(mg.project.discount_rate, T)
    # Generator diagnostics
    Pfc = value.(md["Pfc"])
    Efc = value(md["Efc"])
    power_rated_fc = value(md["power_rated_fc"])
    # operation hours/y
    fc_hours = sum(Pfc .> 1e-6*power_rated_fc)*dt*365/ndays
    fc_hours_lin = relax_gain*Efc/power_rated_fc
    fc_lifetime = mg.dispatchables.fuel_cell[1].lifetime_hours / fc_hours
    fc_lifetime_hlin = mg.dispatchables.fuel_cell[1].lifetime_hours / fc_hours_lin
    # Storage diagnostics
    energy_rated_sto = value(md["energy_rated_sto"])
    E_through_sto = value(md["E_through_sto"])
    sto_cycles = value(E_through_sto/(2*energy_rated_sto)) # c/year
    sto_lifetime_cycles = mg.storage.lifetime_cycles / sto_cycles
    sto_lifetime = min(sto_lifetime_cycles, mg.storage.lifetime_calendar)

    diagnostics = (
        fuel_cell = (
            cost_share = value(md["Cfc"]/md["Cann"]),
            cost_share_fuel = value(md["Cfc_combustible"]/md["Cann"]),
            energy = Efc,
            load_share = Efc/Eload_desired,
            hours = fc_hours,
            hours_lin = fc_hours_lin,
            lifetime = fc_lifetime,
            lifetime_hlin = fc_lifetime_hlin,
            CRF = CRFproj(fc_lifetime),
            CRF_hlin = CRFproj(fc_lifetime_hlin),
            CRF_hlin_pwl = value(md["Pfc_rated_ann"] / power_rated_fc)
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
        fc_hours_assum = 2000.,
        relax_gain = 1.0,
        z_tan = [0.20, 0.28, 0.37, 0.50, 0.68, 1.0, 1.7, 4.0],
        model_custom=nothing,
        H2_price = 0.0,
        infinite_storage::Bool = false,
    )

    # Setup optimization model
    mg, model_data = setup_optim_mg_jump(optimizer;
        shed_max, ndays, fixed_lifetimes, fc_hours_assum, relax_gain, z_tan,
        create_mg_base, H2_price, infinite_storage
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
    Ctank = infinite_storage ? 0.0 : value(md["Ctank"])
    xopt = value.([
        md["power_rated_fc"]
        md["power_rated_ele"]
        md["energy_rated_sto"]
        Ctank
        md["power_rated_pv"]
        md["power_rated_wind"]
    ])

    diagnostics = diagnostics_mg_jump(mg, model_data, ndays, relax_gain)

    K = length(md["Pnl"])
    traj = (
        Pfc = value.(md["Pfc"]),
        Pgen = zeros(K),
        Pele = value.(md["Pele"]),
        Phb   = zeros(K),
        Psto = value.(md["Psto_dis"] - md["Psto_cha"]),
        Esto = value.(md["Esto"]),
        Ctank = Ctank,
        LoH = zeros(K),
        LoF   = zeros(K),
        Pcurt = value.(md["Pspill"]),
        Pdump = zeros(K)
    )
    status = termination_status(model)
    return xopt, LCOE_opt, diagnostics, traj, md, mg, status
end