# Parameters and time series for a Microgrid project with *wind* and *solar* sources,
# plus a *battery* and a *dispatchable generator*.
#
# Values gathered from the Microgrid_Wind-Solar.ipynb notebook,
# to be used in the sizing optimization notebook.
#
# July 2024 update taking advantage of component mutability since v0.11.0

using Microgrids # for capacity_from_wind
using CSV, DataFrames

println("Base Wind-Solar Microgrid project definition functions...")

"""
load a 1 year times series for microgrid project from `Ouessant_data_2016.csv` table:
- load power request
- capacity factor for solar PV production
- capacity factor for wind power production
returned as the `NamedTuple` `(;Pload, cf_pv, cf_wind)`

should be called as `const tseries = load_microgrid_tseries()`
"""
function load_microgrid_tseries()
    ### Input time series
    fname = "data/Ouessant_data_2016.csv"
    println("loading times series from $fname...")
    data = DataFrame(CSV.File("$(@__DIR__)/$fname"))
    
    # Simulation steps
    nsteps = length(data.Load)
    @assert(nsteps == 8760) # 1 year at an hourly timestep
    
    # Split load, solar and wind data:
    Pload = data.Load # kW
    cf_pv =  data.Ppv1k ./ 1000; # convert to kW/kWp
    wind_speed = data.Wind; # m/s
    
    # Calibrate wind speed data against a mast measurement
    ws_gain = 1.059 # ratio of Mast's mean /PVGIS' mean
    wind_speed = ws_gain*wind_speed
    
    # Generate wind power capacity factor time series from wind speed
    S_D52 = pi * (52/2)^2 # rotor swept area m²
    TSP_D52 = 900e3/S_D52 # W/m²
    v_out = 25.0 # m/s
    Cp_D52, α_D52 = 0.521, 3.1; # fitted from actual power curve
    cf_wind = capacity_from_wind.(wind_speed; TSP=TSP_D52, Cp=Cp_D52, v_out=v_out, α=α_D52)

    return (;Pload, cf_pv, cf_wind)
end

"""
Create a base Microgrid project description

The Microgrid project comes with *wind* and *solar* sources,
plus a *battery* and a *dispatchable generator*,
created with some default parameters

The function takes as an input the global variable `tseries`,
with the fields `Pload`, `cf_pv`, `cf_wind`, as created by `load_microgrid_tseries()`.
"""
function create_mg_base()
    ### Microgrid parameters
    # Project parameters
    lifetime = 25 # yr
    discount_rate = 0.05
    timestep = 1. # h
    
    # Parameters chosen to be common to all Components (but they could differ)
    replacement_price_ratio = 1.0
    salvage_price_ratio = 1.0
    
    # Dispatchable generator (Diesel)
    power_rated_gen = 1800.  # rated power (kW) → /2 to see some load shedding
    # fuel parameters:
    fuel_intercept = 0.0 # fuel curve intercept (L/h/kW_max)
    fuel_slope = 0.240 # fuel curve slope (L/h/kW)
    # prices:
    fuel_price = 1. # fuel price ($/L)
    investment_price_gen = 400. # initial investiment price ($/kW)
    om_price_gen = 0.02 # operation & maintenance price ($/kW/h of operation)
    # other technical parameters
    lifetime_gen = 15000. # generator lifetime (h)
    load_ratio_min = 0.0 # minimum load ratio ∈ [0,1]
    fuel_unit = "L" # fuel quantity unit (used in fuel price and consumption curve parameters)
    
    # Battery energy storage
    energy_rated_sto = 5000. # rated energy capacity (kWh)
    investment_price_sto = 350. # initial investiment price ($/kWh)
    om_price_sto = 10. # operation and maintenance price ($/kWh/y)
    lifetime_sto = 15. # calendar lifetime (y)
    lifetime_cycles = 3000. # maximum number of cycles over life (1)
    # Parameters with default values
    charge_rate = 1.0 # max charge power for 1 kWh (kW/kWh = h^-1)
    discharge_rate = 1.0 # max discharge power for 1 kWh (kW/kWh = h^-1)
    loss_factor_sto = 0.05 # linear loss factor α (round-trip efficiency is about 1 − 2α) ∈ [0,1]
    SoC_min = 0. # minimum State of Charge ∈ [0,1]
    SoC_ini = 0. # initial State of Charge ∈ [0,1]
    
    # Photovoltaic (PV) generation
    power_rated_pv = 3000. # rated power (kW)
    irradiance = tseries.cf_pv # global solar irradiance incident on the PV array (kW/m²)
    investment_price_pv = 1200. # initial investiment price ($/kW)
    om_price_pv = 20.# operation and maintenance price ($/kW/y)
    lifetime_pv = 25. # lifetime (y)
    # Parameters with default values
    derating_factor_pv = 1.0 # derating factor (or performance ratio) ∈ [0,1]
    
    # Wind power generation
    power_rated_wind = 900. # rated power (kW)
    cf_wind = tseries.cf_wind
    investment_price_wind = 3500. # initial investiment price ($/kW)
    om_price_wind = 100.# operation and maintenance price ($/kW/y)
    lifetime_wind = 25. # lifetime (y)
    
    ### Microgrid project description structures
    project = Project(lifetime, discount_rate, timestep, "€")
    
    # Create components
    gen = DispatchableGenerator(power_rated_gen,
        fuel_intercept, fuel_slope, fuel_price,
        investment_price_gen, om_price_gen, lifetime_gen,
        load_ratio_min,
        replacement_price_ratio, salvage_price_ratio, fuel_unit)
    batt = Battery(energy_rated_sto,
        investment_price_sto, om_price_sto, lifetime_sto, lifetime_cycles,
        charge_rate, discharge_rate, loss_factor_sto, SoC_min, SoC_ini,
        replacement_price_ratio, salvage_price_ratio)
    pv = Photovoltaic(power_rated_pv, irradiance,
        investment_price_pv, om_price_pv,
        lifetime_pv, derating_factor_pv,
        replacement_price_ratio, salvage_price_ratio)
    windgen = WindPower(power_rated_wind, cf_wind,
        investment_price_wind, om_price_wind,
        lifetime_wind,
        replacement_price_ratio, salvage_price_ratio)
    
    mg = Microgrid(project, tseries.Pload, gen, batt, [pv, windgen])

    return mg
end
