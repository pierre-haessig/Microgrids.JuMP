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

    #Fuel cell
    power_rated_fc = 1800. # rated power capacity (KW)
    cons_intercept_fc=0. #
    cons_rate_fc = 0.0625 # consumption rate (KgH2/KWhé)
    cons_price_fc = 0. #
    investment_price_fc = 1600. # initial investment price  ($/KW)
    om_price_fc = 45. # operation and maintenance price ($/kW/y)
    om_price_hour_fc = 0. # operation and maintenance price ($/kW/h)
    lifetime_fc_y = 15.
    lifetime_fc_h = 45000. # Fuel Cell lifetime (h)
    lifetime_fc_starts = 5000. # fuel cell maximum starts on
    load_min_ratio_fc = 0.05 # minimum load ratio ∈ [0,1]
    input_unit_fc= "Kg"
    output_unit_fc="KWh"

    #hy tank
    capacity_rated_hytank = 12000. # rated power capacity (kg)
    investment_price_hytank = 400. # initial investment price  ($/kg)
    hy_price = 14. # initial hydrogen price ($/kg)
    om_price_hytank = 4. # operation and maintenance price ($/kg/y)
    lifetime_hytank = 25. # calendar lifetime (y)
    loss_factor_hytank = 0. # hydrogen used on site 
    LoH_ini_ratio = 0.0 # Initial load ratio ∈ [0,1]
    LoH_min_ratio = 0.0 # minimum load ratio ∈ [0,1]
    LoH_max_ratio = 1. # maximum load ratio ∈ [0,1]

    #Haber Bosch numbers to update!!! they have been taken arbitrary
    power_rated_hb = 4000. # rated power capacity (KW)
    cons_intercept_hb=0. #
    cons_rate_hb = 0.344 # consumption rate (KWhé/KgH2)
    cons_price_hb = 0. #
    investment_price_hb_only = 0.418*580 #energy consumption of HB process for 1 kW of HB+ASU multiplied by CAPEX of HB ($/kW) [Achour et al., 2025]
    investment_price_asu_only = 0.582*224 #energy consumption of ASU for 1 kW of HB+ASU multiplied by CAPEX of ASU ($/kW) [Achour et al., 2025]
    investment_price_hb = investment_price_asu_only+investment_price_hb_only # initial investment price ($/KW)
    om_price_hb = 0.017*1000 # operation and maintenance price ($/kW/y)[Achour et al., 2025]
    om_price_hour_hb = 0. # operation and maintenance price ($/kW/h)
    lifetime_hb_y = 25.
    lifetime_hb_h = 45000. # FOR NOW NOT USED. THIS IS THE VALUE FOR FUELL CELL.
    lifetime_hb_starts = 5000. # fuel cell maximum starts on
    load_min_ratio_hb= 0.3 # minimum load ratio ∈ [0,1]
    input_unit_hb= "Kwh"
    output_unit_hb="Kg"



    # Battery energy storage
    energy_rated_sto = 5000. # rated energy capacity (kWh)
    investment_price_sto = 350. # initial investiment price ($/kWh)
    om_price_sto =10. # operation and maintenance price ($/kWh/y)
    lifetime_sto = 15. # calendar lifetime (y)
    lifetime_cycles = 3000. # maximum number of cycles over life (1)
    # Parameters with default values
    charge_rate = 1.0 # max charge power for 1 kWh (kW/kWh = h^-1)
    discharge_rate = 1.0 # max discharge power for 1 kWh (kW/kWh = h^-1)
    loss_factor_sto = 0.05 # linear loss factor α (round-trip efficiency is about 1 − 2α) ∈ [0,1]
    SoC_min = 0. # minimum State of Charge ∈ [0,1]
    SoC_ini = 0. # initial State of Charge ∈ [0,1]
    SoC_max = 1. # initial State of Charge ∈ [0,1]

    #electrolyzer
    power_rated_elyz = 2000. # rated power capacity (kW)
    cons_intercept_elyz= 0. # consumption rate (KWhé/KgH2)
    cons_rate_elyz = 56.  # consumption rate (KWhé/KgH2)
    cons_price_elyz = 0.#
    investment_price_elyz = 1600. # initial investment price  ($/kW)
    om_price_elyz = 44. # operation and maintenance price ($/kW/y)
    om_price_hour_elyz = 0. # operation and maintenance price ($/kW/h)
    lifetime_elyz_y = 20.
    lifetime_elyz_h = 45000. #Electrolyzer lifetime (h)
    lifetime_elyz_starts = 5000. #Electrolyzer maximum starts on 
    load_min_ratio_elyz = 0.05  # minimum load ratio ∈ [0,1]

    input_unit_elyz= "KWh"
    output_unit_elyz="Kg"

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

    #!!!!!!!Haber Bosch with 0 power rated to satisfy microgrid struct!!!
    power_rated_hb = 0. # rated power capacity (KW)
    cons_intercept_hb=0. #
    cons_rate_hb = 0.344 # consumption rate (KWhé/KgH2)
    cons_price_hb = 0. #
    investment_price_hb_only = 0.418*580 #energy consumption of HB process for 1 kW of HB+ASU multiplied by CAPEX of HB ($/kW) [Achour et al., 2025]
    investment_price_asu_only = 0.582*224 #energy consumption of ASU for 1 kW of HB+ASU multiplied by CAPEX of ASU ($/kW) [Achour et al., 2025]
    investment_price_hb = investment_price_asu_only+investment_price_hb_only # initial investment price ($/KW)
    om_price_hb = 0.017*1000 # operation and maintenance price ($/kW/y)[Achour et al., 2025]
    om_price_hour_hb = 0. # operation and maintenance price ($/kW/h)
    lifetime_hb_y = 25.
    lifetime_hb_h = 45000. # FOR NOW NOT USED. THIS IS THE VALUE FOR FUELL CELL.
    lifetime_hb_starts = 5000. # fuel cell maximum starts on
    load_min_ratio_hb= 0.3 # minimum load ratio ∈ [0,1]
    input_unit_hb= "Kwh"
    output_unit_hb="Kg"
    
    ### Microgrid project description structures
    project = Project(lifetime, discount_rate, timestep, "€")
    
    fuel_cell = ProductionUnit(
      power_rated_fc,
      cons_intercept_fc, 
      cons_rate_fc,
      cons_price_fc,
      investment_price_fc, 
      om_price_hour_fc, 
      om_price_fc,
      lifetime_fc_y,
      lifetime_fc_h,
      lifetime_fc_starts,
      load_min_ratio_fc,
      replacement_price_ratio, 
      salvage_price_ratio,
      input_unit_fc,
      output_unit_fc
      )
    electrolyzer = ProductionUnit(
      power_rated_elyz,
      cons_intercept_elyz,
      cons_rate_elyz,
      cons_price_elyz, 
      investment_price_elyz, 
      om_price_hour_elyz,
      om_price_elyz, 
      lifetime_elyz_y,
      lifetime_elyz_h,
      lifetime_elyz_starts,
      load_min_ratio_elyz,
      replacement_price_ratio,
      salvage_price_ratio,
      input_unit_elyz,
      output_unit_elyz
    )
    batt = Battery(
      energy_rated_sto,
      investment_price_sto,
      om_price_sto, 
      lifetime_sto, 
      lifetime_cycles,
      charge_rate, 
      discharge_rate, 
      loss_factor_sto, 
      SoC_min, 
      SoC_max,
      SoC_ini,
      replacement_price_ratio, 
      salvage_price_ratio
    )
    pv = Photovoltaic(
      power_rated_pv, 
      irradiance,
      investment_price_pv, 
      om_price_pv,
      lifetime_pv, 
      derating_factor_pv,
      replacement_price_ratio, 
      salvage_price_ratio
    )
    windgen = WindPower(
      power_rated_wind, 
      cf_wind,
      investment_price_wind, 
      om_price_wind,
      lifetime_wind,
      replacement_price_ratio, 
      salvage_price_ratio
    )

    haber = ProductionUnit(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "", "")#no haber-bosh in our case
    tank_hy = Tank(capacity_rated_hytank, investment_price_hytank, om_price_hytank,lifetime_hytank, loss_factor_hytank,LoH_ini_ratio, LoH_min_ratio, LoH_max_ratio, hy_price, 0.0, 0.0)
    tank_0 = Tank(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0) #infinite reservoir, so we don't need them
    gen = ProductionUnit(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "", "")#no generaator in our model
    mg = Microgrid(
    project,
    tseries.Pload,
    DispatchableCompound([gen],[fuel_cell]),
    [electrolyzer],
    haber,
    TankCompound(tank_hy, tank_0),
    batt,
    [pv, windgen])

    return mg
end