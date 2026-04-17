using Microgrids
using JuMP
using HiGHS
using DataFrames
using CSV

include("../data/Microgrid_model_data.jl")
include("economics.jl")
include("simulate.jl")
include("create_microgrid.jl")

import MathOptInterface as MOI

const tseries = load_microgrid_tseries()

function LCOE_model(H2_price, td, infinite_storage::Bool = false)

    xopt, LCOE_opt, diagnostics, traj, md, mg, status =
        optim_mg_jump(optimizer, H2_price = H2_price, create_mg_base = create_mg_base, infinite_storage = infinite_storage)


    if infinite_storage
        H2_per_kWh = 1 / mg.electrolyzer[1].consumption_slope
        H2_ele = value(md["Eele"]) * H2_per_kWh
        H2_fc = value(md["Efc"]) * mg.dispatchables.fuel_cell[1].consumption_slope

        Q, range_Q, cost_Q = Q_hydro_overtime(mg, md, td)

        C_combustible_net = H2_price * H2_fc - H2_price * H2_ele

        real_Cann = value(md["Cann"]) - C_combustible_net + cost_Q
        real_LCOE = real_Cann / value(md["Eload_desired"])

        return xopt, status, Q, cost_Q / value(md["Eload_desired"]), LCOE_opt, real_LCOE
    end

    return xopt, status, LCOE_opt
end

function Study(H2_price_range = [], infinite_storage::Bool = false)

    LCOE_array = Float64[]
    x_array = Vector{Any}()
    status_array = Vector{MOI.TerminationStatusCode}()

    if infinite_storage
        real_LCOE_array = Float64[]
        Q_array = Vector{Vector{Float64}}()
        cost_Q_array = Float64[]

        K = 365 * 24
        td = collect((0:K-1) / 24)

        for i in H2_price_range
            xopt, status, Q, cost_Q, LCOE, real_LCOE =
                LCOE_model(i, td, true)

            push!(LCOE_array, LCOE)
            push!(real_LCOE_array, real_LCOE)
            push!(x_array, xopt)
            push!(status_array, status)
            push!(Q_array, Q)
            push!(cost_Q_array, cost_Q)
        end

        return x_array, status_array, Q_array, cost_Q_array, real_LCOE_array, LCOE_array
    end

    # single run (H2_price irrelevant)
    xopt, status, LCOE = LCOE_model(0, [], false)

    push!(LCOE_array, LCOE)
    push!(x_array, xopt)
    push!(status_array, status)

    return x_array, status_array, LCOE_array
end

function dataframe_result(result, H2_price_range, infinite_storage::Bool = false)

    if infinite_storage
        x_array, status_array, Q_array, cost_Q_array, real_LCOE_array, LCOE_array = result
        df_result = DataFrame(
        H2_price = H2_price_range,
        LCOE = 1000 .* LCOE_array,
        real_LCOE = 1000 .* real_LCOE_array,
        power_rated_fc = getindex.(x_array, 1),
        power_rated_ele = getindex.(x_array, 2),
        energy_rated_sto = getindex.(x_array, 3),
        power_rated_pv = getindex.(x_array, 5),
        power_rated_wind = getindex.(x_array, 6),
        status = string.(status_array),
        cost_Q = 1000 .* cost_Q_array
    )

        df_Q = DataFrame(
            H2_price = Float64[],
            time = Int[],
            Q = Float64[]
        )

        for i in eachindex(H2_price_range)
            Q_i = Q_array[i]
            for t in eachindex(Q_i)
                push!(df_Q, (
                    H2_price = H2_price_range[i],
                    time = t,
                    Q = Q_i[t]
                ))
            end
        end

        return df_result, df_Q
    end

    x_array, status_array, LCOE_array = result

    df_result = DataFrame(
        LCOE = 1000 .* LCOE_array,
        power_rated_fc = getindex.(x_array, 1),
        power_rated_ele = getindex.(x_array, 2),
        energy_rated_sto = getindex.(x_array, 3),
        LoH_rated = getindex.(x_array, 4),
        power_rated_pv = getindex.(x_array, 5),
        power_rated_wind = getindex.(x_array, 6),
        status = string.(status_array),
    )
    return df_result
    
end

function export_to_csv(df::DataFrame; prefix="result_for_H2_price")

    if :H2_price in names(df)
        h2_min = minimum(df.H2_price)
        h2_max = maximum(df.H2_price)
        filename = "plots/CSV_results/$(prefix)_$(h2_min)_to_$(h2_max).csv"
        CSV.write(filename, df)
        return filename
    end
    
    filename = "plots/CSV_results/no_infinite storage"
    CSV.write(filename, df)
    return filename
end