using Microgrids
using JuMP
using HiGHS
using DataFrames
using CSV

include("example/Microgrid_model_data.jl")
include("Microgrid_JuMP_common.jl")

function LCOE_bf_af(H2_price,td)
    xopt, LCOE_opt, diagnostics, traj, md, mg =
            optim_mg_jump(optimizer, H2_price = H2_price, create_mg_base = create_mg_base)

        H2_per_kWh = 1 / mg.electrolyzer[1].consumption_slope
        H2_ele = value(md["Eele"]) * H2_per_kWh
        H2_fc = value(md["Efc"]) * mg.dispatchables.fuel_cell[1].consumption_slope

        Q, range_Q, cost_Q = Q_hydro_overtime(mg, md, td)
        C_combustible_net = H2_price * H2_fc - H2_price * H2_ele
        real_Cann = value(md["Cann"]) 
                     - C_combustible_net
                     + cost_Q

        LCOE = 1000 * value(md["Cann"]) / value(md["Eload_desired"])
        real_LCOE = 1000 * real_Cann / value(md["Eload_desired"])
    
        println("Eele = ", value(md["Eele"]))
        println("Efc  = ", value(md["Efc"]))

        println("H2_ele = ", H2_ele)
        println("H2_fc  = ", H2_fc)

        println("ΔH2 = ", H2_fc - H2_ele)

        println("cost_Q = ", cost_Q)
        println("C_combustible_net = ", C_combustible_net)
    return xopt, LCOE, real_LCOE
end

function Study_H2_price(H2_price_range)
    LCOE_H2_price = Float64[]
    LCOE_bf = Float64[]

    K = 365*24
    td = collect((0:K-1)/24)

    for i in H2_price_range
        xopt, LCOE, real_LCOE = LCOE_bf_af(i, td)
        push!(LCOE_H2_price, real_LCOE)
        push!(LCOE_bf, LCOE)
    end

    return LCOE_H2_price, LCOE_bf
end


function dataframe_H2(LCOE_tuple, H2_price_range)
    real_LCOE_array, LCOE_array = LCOE_tuple

    df_H2 = DataFrame(
        H2_price = Float64[],
        LCOE = Float64[],
        real_LCOE = Float64[]
    )

    for i in eachindex(H2_price_range)
        push!(df_H2, (
            H2_price = H2_price_range[i],
            LCOE = real_LCOE_array[i],
            real_LCOE = LCOE_array[i]
        ))
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


"H2_price_range = 0:1:30
results = Study_H2_price(H2_price_range)
df = dataframe_H2(results, H2_price_range)
export_to_csv(df)"
K = 365*24
td = collect((0:K-1)/24)
x, L, real_L = LCOE_bf_af(14, td)
println(L)
println(real_L)
println(x)