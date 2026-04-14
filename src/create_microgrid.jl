optimizer = HiGHS.Optimizer
"""mg = create_mg_base()"""

"""Pload_max = maximum(tseries.Pload) # kW
Pload = tseries.Pload
power_rated_gen_max = 0
power_rated_fc_max = 1.2 * Pload_max
power_rated_ele_max = 1.2* Pload_max
power_rated_hb_max = 0
energy_rated_sto_max = 10.0 * Pload_max
capacity_max = 0
power_rated_pv_max = 1.0 * Pload_max
power_rated_wind_max = 5.0 * Pload_max"""


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
    tank_hy = Tank(
        mg.tanks.h2Tank.capacity, 
        mg.tanks.h2Tank.investment_price, 
        mg.tanks.h2Tank.om_price,
        mg.tanks.h2Tank.lifetime, 
        mg.tanks.h2Tank.loss_factor,
        mg.tanks.h2Tank.ini_filling_ratio, 
        mg.tanks.h2Tank.min_filling_ratio, 
        mg.tanks.h2Tank.max_filling_ratio, 
        mg.tanks.h2Tank.combustible_price, 
        0.0, 
        0.0)

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
        TankCompound(tank_hy, tank_hy),
        batt_new,
        [pv_new, windgen_new]
    )

    return mg_new
end