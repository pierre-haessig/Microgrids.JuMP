# Common code for Microgrid definition, simulation and optimization
# Pierre Haessig, September-October 2024
# This is a script and function-wrapped version of the code developped in Microgrid_optimization_JuMP.ipynb
# which is then reused in other notebooks

using Microgrids
using JuMP
using HiGHS
using DataFrames
include("example/Microgrid_model_data.jl")
const tseries = load_microgrid_tseries()

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
println("- Q_hydro_overtime")

optimizer = HiGHS.Optimizer
mg = create_mg_base()

Pload_max = maximum(tseries.Pload) # kW
Pload = tseries.Pload
power_rated_gen_max = 0
power_rated_fc_max = 1.2 * Pload_max
power_rated_ele_max = 1.2* Pload_max
power_rated_hb_max = 0
energy_rated_sto_max = 10.0 * Pload_max
capacity_max = 0
power_rated_pv_max = 1.0 * Pload_max
power_rated_wind_max = 5.0 * Pload_max


function create_microgrid(x; create_mg_base=create_mg_base)
    mg = create_mg_base()

    irr_raw = mg.nondispatchables[1].irradiance
    irr_norm = irr_raw ./ maximum(irr_raw)

    cf_raw = mg.nondispatchables[2].capacity_factor
    cf_norm = cf_raw ./ maximum(cf_raw)

    gen_new = ProductionUnit(0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, Inf, Inf, Inf, 1.0, 1.0, 1.0, "", "")

    fuel_cell_new = ProductionUnit(x[1],
        mg.dispatchables.fuel_cell[1].consumption_intercept,
        mg.dispatchables.fuel_cell[1].consumption_slope,
        mg.dispatchables.fuel_cell[1].combustible_price,
        mg.dispatchables.fuel_cell[1].investment_price,
        mg.dispatchables.fuel_cell[1].om_price_hourly,
        mg.dispatchables.fuel_cell[1].om_price,
        mg.dispatchables.fuel_cell[1].lifetime_calendar,
        mg.dispatchables.fuel_cell[1].lifetime_hours,
        mg.dispatchables.fuel_cell[1].lifetime_on_off,
        mg.dispatchables.fuel_cell[1].minimum_load_ratio,
        mg.dispatchables.fuel_cell[1].replacement_price_ratio,
        mg.dispatchables.fuel_cell[1].salvage_price_ratio,
        mg.dispatchables.fuel_cell[1].input_unit,
        mg.dispatchables.fuel_cell[1].output_unit)

    electrolyzer_new = ProductionUnit(x[2],
        mg.electrolyzer[1].consumption_intercept,
        mg.electrolyzer[1].consumption_slope,
        mg.electrolyzer[1].combustible_price,
        mg.electrolyzer[1].investment_price,
        mg.electrolyzer[1].om_price_hourly,
        mg.electrolyzer[1].om_price,
        mg.electrolyzer[1].lifetime_calendar,
        mg.electrolyzer[1].lifetime_hours,
        mg.electrolyzer[1].lifetime_on_off,
        mg.electrolyzer[1].minimum_load_ratio,
        mg.electrolyzer[1].replacement_price_ratio,
        mg.electrolyzer[1].salvage_price_ratio,
        mg.electrolyzer[1].input_unit,
        mg.electrolyzer[1].output_unit)

    haber = ProductionUnit(0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, Inf, Inf, Inf, 1.0, 1.0, 1.0, "", "")

    batt_new = Battery(x[3],
        mg.storage.investment_price,
        mg.storage.om_price,
        mg.storage.lifetime_calendar,
        mg.storage.lifetime_cycles,
        mg.storage.charge_rate,
        mg.storage.discharge_rate,
        mg.storage.loss_factor,
        mg.storage.SoC_min,
        mg.storage.SoC_max,
        mg.storage.SoC_ini,
        mg.storage.replacement_price_ratio,
        mg.storage.salvage_price_ratio)

    tank_0 = Tank(0.0, 1.0, 1.0, Inf, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

    pv_new = Photovoltaic(x[4],
        irr_norm,
        mg.nondispatchables[1].investment_price,
        mg.nondispatchables[1].om_price,
        mg.nondispatchables[1].lifetime,
        mg.nondispatchables[1].derating_factor,
        mg.nondispatchables[1].replacement_price_ratio,
        mg.nondispatchables[1].salvage_price_ratio)

    windgen_new = WindPower(x[5],
        cf_norm,
        mg.nondispatchables[2].investment_price,
        mg.nondispatchables[2].om_price,
        mg.nondispatchables[2].lifetime,
        mg.nondispatchables[2].replacement_price_ratio,
        mg.nondispatchables[2].salvage_price_ratio)

    mg_new = Microgrid(
        mg.project,
        tseries.Pload,
        DispatchableCompound([fuel_cell_new, gen_new], [fuel_cell_new, gen_new]),
        [electrolyzer_new],
        haber,
        TankCompound(tank_0, tank_0),
        batt_new,
        [pv_new, windgen_new]
    )

    return mg_new
end

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
    `shed_max`, `ndays`, `fixed_lifetimes`, `fc_hours_assum`, `relax_gain` and `z_tan`

⚠ Sizing maximum bounds are taken as global variables:
- `power_rated_fc_max`
- `power_rated_ele`
- `energy_rated_sto_max`
- `power_rated_pv_max`
- `power_rated_wind_max`

Also, a global `ts_reduction(x, ndays)` function is needed.
"""
function build_optim_mg_stage!(mg::Microgrid, model_data::Dict{String,Any}, H2_price = 0.0)
    # Retrieve model parameters from the Dict:
    md = model_data # shortcut Dict name
    shed_max = md["shed_max"]
    ndays = md["ndays"]
    fixed_lifetimes = md["fixed_lifetimes"]
    fc_hours_assum = md["fc_hours_assum"]
    relax_gain = md["relax_gain"]
    z_tan = md["z_tan"]

    #println("Building stage problem with $ndays days...")
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
    
    raw_pv = production(mg.nondispatchables[1])
    cf_pv = (raw_pv ./ maximum(raw_pv)) |> ts_reduction_ndays

    raw_wind = production(mg.nondispatchables[2])
    cf_wind = (raw_wind ./ maximum(raw_wind)) |> ts_reduction_ndays

    ### JuMP model definition
    model = md["model"] # JuMP model


    ##  Sizing variables
    md["power_rated_fc"]  = @variable(model, 0 <= power_rated_fc <= power_rated_fc_max)
    md["power_rated_ele"] = @variable(model, 0 <= power_rated_ele <= power_rated_ele_max)
    md["energy_rated_sto"] = @variable(model, 0 <= energy_rated_sto <= energy_rated_sto_max)
    md["power_rated_pv"]   = @variable(model, 0 <= power_rated_pv <= power_rated_pv_max)
    md["power_rated_wind"] = @variable(model, 0 <= power_rated_wind <= power_rated_wind_max)



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

    # Fuel_cells
    md["Pfc"] = @variable(model, Pfc[1:K] >= 0)
    @constraint(model, Pfc .<= power_rated_fc)

    #Electrolyzer
    md["Pele"] = @variable(model, Pele[1:K] >= 0)
    @constraint(model, Pele .<= power_rated_ele)
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
    @constraint(model, Psto_dis .<= mg.storage.discharge_rate * energy_rated_sto) # double the solving time with HiGHS and MOI 1.30.0, when there is the Pfc penalty!!
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
        Pfc - Pele + (Psto_dis - Psto_cha) - Pspill .== Pnl - Pshed,
    )

    if shed_max == 0.0
        #println("zero load shedding allowed")
        fix.(Pshed, 0.0; force=true);
    else
        #println("load shedding allowed: $shed_max")
        Eshed = sum(Pshed)*dt * 365/ndays
        @constraint(model, Eshed <= Eload_desired*shed_max)
    end

    ### Costs
    # Fuel cell costs
    md["Pfc_rated_ann"] = @variable(model, Pfc_rated_ann >= 0) # annualized size
    md["Efc"] = Efc = sum(Pfc)*dt * 365/ndays # Generator yearly energy

    if fixed_lifetimes
        #print("Fixed fuel_cell lifetime hypothesis: ")
        fc_lifetime = mg.dispatchables.fuel_cell[1].lifetime_hours / fc_hours_assum # years
        #println("$fc_lifetime y, assuming $fc_hours_assum  h/y of usage")
        md["cons_Pfc_rated_ann_CRFfc_lifetime"] = @constraint(model,
            Pfc_rated_ann == power_rated_fc * CRFproj(fc_lifetime))
    else
        #println("Usage-dependent fuel_cell lifetime model (relax_gain=", relax_gain,")")
        md["Ufc"] = @variable(model, Ufc >= 0) # cumulated usage
        #println("Xann constraints with z_tan=$z_tan")
        @constraint(model, Ufc == Efc*relax_gain/mg.dispatchables.fuel_cell[1].lifetime_hours); # kW/y
        cpwl_fc = cons_Xann_usage!(model,
            Pfc_rated_ann, power_rated_fc, Ufc,
            discount_rate, z_tan)
    end
    # question for O&M wit fixed lifetime: use assumed fixed fc hours ?
    # or use the convexified hours (Efc*relax_gain) like for Ufc?
    md["Cfc_om"] = Cfc_om = fixed_lifetimes ?
        mg.dispatchables.fuel_cell[1].om_price_hourly * fc_hours_assum * power_rated_fc :
        mg.dispatchables.fuel_cell[1].om_price_hourly * relax_gain * Efc

    md["Cfc_combustible"] = Cfc_combustible = H2_price * mg.dispatchables.fuel_cell[1].consumption_slope * Efc;# $/y
    md["Cfc"] = Cfc = mg.dispatchables.fuel_cell[1].investment_price * Pfc_rated_ann +
                        Cfc_om + 
                        Cfc_combustible # $/y

    # Electroolyzer costs
    md["Pele_rated_ann"] = @variable(model, Pele_rated_ann >= 0) # annualized size
    md["Eele"] = Eele = sum(Pele)*dt * 365/ndays # Generator yearly energy

    if fixed_lifetimes
        #print("Fixed electrolyzer lifetime hypothesis: ")
        ele_lifetime = mg.electrolyzer[1].lifetime_hours / ele_hours_assum # years
        #println("$ele_lifetime y, assuming $ele_hours_assum  h/y of usage")
        md["cons_Pele_rated_ann_CRFele_lifetime"] = @constraint(model,
            Pele_rated_ann == power_rated_ele * CRFproj(ele_lifetime))
    else
        #println("Usage-dependent electrolyzer lifetime model (relax_gain=", relax_gain,")")
        md["Uele"] = @variable(model, Uele >= 0) # cumulated usage
        #println("Xann constraints with z_tan=$z_tan")
        @constraint(model, Uele == Eele*relax_gain/mg.electrolyzer[1].lifetime_hours); # kW/y
        cpwl_ele = cons_Xann_usage!(model,
            Pele_rated_ann, power_rated_ele, Uele,
            discount_rate, z_tan)
    end
    # question for O&M wit fixed lifetime: use assumed fixed ele hours ?
    # or use the convexified hours (Eele*relax_gain) like for Uele?
    md["Cele_om"] = Cele_om = fixed_lifetimes ?
        mg.electrolyzer[1].om_price_hourly * ele_hours_assum * power_rated_ele :
        mg.electrolyzer[1].om_price_hourly * relax_gain * Eele

    md["Cele_combustible"] = Cele_combustible = H2_price / mg.electrolyzer[1].consumption_slope * Eele;# $/y
    md["Cele"] = Cele = mg.electrolyzer[1].investment_price * Pele_rated_ann + 
                        Cele_om - 
                        Cele_combustible # $/y
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
    md["Cann"] = Cann = Cfc + Csto +  Cpv + Cwind + Cele
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

    build_optim_mg_stage!(mg, model_data, H2_price)

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
        H2_price = 0.0
    )

    # Setup optimization model
    mg, model_data = setup_optim_mg_jump(optimizer;
        shed_max, ndays, fixed_lifetimes, fc_hours_assum, relax_gain, z_tan,
        create_mg_base, H2_price
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
        md["power_rated_fc"]
        md["power_rated_ele"]
        md["energy_rated_sto"]
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
        LoH   = zeros(K),
        LoF   = zeros(K),
        Pcurt = value.(md["Pspill"]),
        Pdump = zeros(K)
    )
    return xopt, LCOE_opt, diagnostics, traj, md, mg
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

function Q_hydro_overtime(mg::Microgrid, md ,investment_price_hytank::Real=0.0, td::Vector{Float64} = zeros(365*24))
    #Variables needed
    Pele = value.(md["Pele"])
    cons_rate_elyz = mg.electrolyzer[1].consumption_slope
    Pfc   = value.(md["Pfc"])
    cons_rate_fc = mg.dispatchables.fuel_cell[1].consumption_slope
    dt = td[2] - td[1]
    K = length(td)

    Q = cumsum(((Pele[1:K] ./ cons_rate_elyz) .- (Pfc[1:K] .* cons_rate_fc)) .* dt)
    range_Q = maximum(Q) - minimum(Q)
    cost_Q = range_Q * investment_price_hytank * CRF(mg.project.discount_rate, mg.storage.lifetime_calendar)
    return (Q, range_Q, cost_Q)
end

