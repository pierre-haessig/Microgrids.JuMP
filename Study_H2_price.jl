using Microgrids
using JuMP
using HiGHS
using DataFrames
using CSV

include("example/Microgrid_model_data.jl")
include("Microgrid_JuMP_common.jl")

function Study_H2_price(H2_price_range)
    LCOE_H2_price = []
    K = 365*24
    td = collect((0:K-1)/24)

    for i in H2_price_range

        xopt, LCOE_opt, diagnostics, traj, md, mg =
            optim_mg_jump(optimizer, H2_price = i, create_mg_base = create_mg_base)

        H2_per_kWh = 1 / mg.electrolyzer[1].consumption_slope
        H2_prod = value(md["Eele"]) * H2_per_kWh
        H2_fc = value(md["Efc"]) * mg.dispatchables.fuel_cell[1].consumption_slope

        Q, range_Q, cost_Q = Q_hydro_overtime(mg, md, 400, td)

        real_Cann =  value(md["Cann"]) +
                     i * H2_fc -
                     i * H2_prod +
                     cost_Q

        real_LCOE = 1000*real_Cann / value(md["Eload_desired"])
        push!(LCOE_H2_price, real_LCOE)
    end

    return LCOE_H2_price
end

function dataframe_H2(LCOE_H2_price, H2_price_range)
    df_H2 = DataFrame(
        H2_price = Float64[],
        LCOE = Float64[]
        
    )

    for i in eachindex(LCOE_H2_price)
        push!(df_H2, (H2_price = H2_price_range[i], LCOE = LCOE_H2_price[i]))
    
    end

    return df_H2
end

function export_to_csv(df::DataFrame; prefix="LCOE_for_H2_price")
    h2_min = minimum(df.H2_price)
    h2_max = maximum(df.H2_price)
    filename = "$(prefix)_$(h2_min)_to_$(h2_max).csv"

    CSV.write(filename, df)
    return nothing
end


H2_price_range = range(0, 15; step = 0.1)
export_to_csv(dataframe_H2(Study_H2_price(H2_price_range), H2_price_range))