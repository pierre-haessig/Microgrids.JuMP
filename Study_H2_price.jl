using Microgrids
using JuMP
using HiGHS
using DataFrames

include("example/Microgrid_model_data.jl")
include("Microgrid_JuMP_common.jl")
const tseries = load_microgrid_tseries()


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
        print(LCOE_H2_price)
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
 
function export_to_csv(df::DataFrame, filename::String="result.csv")
    CSV.write(filename, df)
    return nothing
end


H2_price_range = range(100, 1000; step = 100)
export_to_csv(dataframe_H2(Study_H2_price(H2_price_range), H2_price_range))