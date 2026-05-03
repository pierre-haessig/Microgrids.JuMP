using Microgrids
using JuMP
using HiGHS
using DataFrames
using CSV
using Base.Threads

Threads.nthreads() #print number of threads used

include("../data/Microgrid_model_data.jl")
include("economics.jl")
include("simulate.jl")
include("create_microgrid.jl")

import MathOptInterface as MOI

const tseries = load_microgrid_tseries()

function LCOE_model(H2_price, td, infinite_storage::Bool = false)

    xopt, LCOE_opt, diagnostics, traj, md_opt, mg_opt, status =
        optim_mg_jump(optimizer, H2_price = H2_price, create_mg_base = create_mg_base, infinite_storage = infinite_storage)


    if infinite_storage
        H2_per_kWh = 1 / mg_opt.electrolyzer[1].consumption_slope
        H2_ele = value(md_opt["Eele"]) * H2_per_kWh
        H2_fc = value(md_opt["Efc"]) * mg_opt.dispatchables.fuel_cell[1].consumption_slope

        Q, range_Q, cost_Q = Q_hydro_overtime(mg_opt, md_opt, td)

        C_combustible_net = H2_price * H2_fc - H2_price * H2_ele

        real_Cann = value(md_opt["Cann"]) - C_combustible_net + cost_Q
        real_LCOE = real_Cann / value(md_opt["Eload_desired"])

        return xopt, status, Q, cost_Q / value(md_opt["Eload_desired"]), range_Q, LCOE_opt, real_LCOE, md_opt, mg_opt
    end

    return xopt, status, LCOE_opt, md_opt, mg_opt
end

function Study(H2_price_range = [], infinite_storage::Bool = false)

    N = length(H2_price_range)

    LCOE_array = Vector{Float64}(undef, N)
    x_array = Vector{Any}(undef, N)
    status_array = Vector{MOI.TerminationStatusCode}(undef, N)
    md_array = Vector{Any}(undef, N)
    mg_array = Vector{Any}(undef, N)

    if infinite_storage
        real_LCOE_array = Vector{Float64}(undef, N)
        Q_array = Vector{Vector{Float64}}(undef, N)
        cost_Q_array = Vector{Float64}(undef, N)
        range_Q_array = Vector{Float64}(undef, N)

        K = 365 * 24
        td = collect((0:K-1) / 24)

        Threads.@threads for idx in eachindex(H2_price_range)
            i = H2_price_range[idx]

            xopt, status, Q, cost_Q, range_Q, LCOE, real_LCOE, md_opt, mg_opt =
                LCOE_model(i, td, true)

            LCOE_array[idx] = LCOE
            real_LCOE_array[idx] = real_LCOE
            x_array[idx] = xopt
            status_array[idx] = status
            Q_array[idx] = Q
            cost_Q_array[idx] = cost_Q
            range_Q_array[idx] = range_Q
            md_array[idx] = md_opt
            mg_array[idx] = mg_opt
        end

        return x_array, status_array, Q_array, cost_Q_array,
               range_Q_array, real_LCOE_array, LCOE_array,
               md_array, mg_array
    end
end

function dataframe_result(result, H2_price_range, infinite_storage::Bool = false)

    if infinite_storage
        x_array, status_array, Q_array, cost_Q_array, range_Q_array, real_LCOE_array, LCOE_array, md_array, mg_array = result
        @show fieldnames(typeof(mg_array[1]))
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
        cost_Q = 1000 .* cost_Q_array,
        range_Q = range_Q_array,
        md = md_array,
        mg = mg_array
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

    x_array, status_array, LCOE_array, md_array, mg_array = result

    df_result = DataFrame(
        LCOE = 1000 .* LCOE_array,
        power_rated_fc = getindex.(x_array, 1),
        power_rated_ele = getindex.(x_array, 2),
        energy_rated_sto = getindex.(x_array, 3),
        LoH_rated = getindex.(x_array, 4),
        power_rated_pv = getindex.(x_array, 5),
        power_rated_wind = getindex.(x_array, 6),
        status = string.(status_array),
        md = md_array,
        mg = mg_array
    )
    return df_result
    
end

function export_to_csv(df::DataFrame; prefix="result_for_H2_price")
    df_clean = select(df, Not(["md", "mg"]))
    if "H2_price" ∈ names(df_clean)
        h2_min = minimum(df_clean.H2_price)
        h2_max = maximum(df_clean.H2_price)
        filename = "../plots/CSV_results/$(prefix)_$(h2_min)_to_$(h2_max).csv"
    else
        filename = "../plots/CSV_results/no_infinite_storage.csv"
    end

    CSV.write(filename, df_clean, bufsize=10_000_000)
    return filename
end