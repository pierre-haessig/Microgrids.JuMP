# ROADMAP


## First file to open

Notebook jupyter that synthesize the results and the overall model :

    "src/H2_optimization.ipynb"

Presentation of the research project to understand the motive and the results we get with this fork :

    "Research_project_presentation.pdf"

To this extent to read it directly in VS code i advise to install the extension "VS code pdf"

## Code outline 



*version april 2026*



### create_microgrid.jl

Functions to create the object microgrid based on the data in :

    "data/Microgrid_model_data.jl"

### economics.jl

Set of functions to calculate the cost of a microgrid biased by artificial and unbiased price

Function to know the actual cost of the H2 tank

    "function Q_hydro_overtime(mg::Microgrid, md, td::Vector{Float64} = zeros(365*24))"

### mg_jump.jl

Set of functions to optimize a microgrid based on the load during a set period of time.

### simulate.jl

Set of functions to simulate the trajectory of the optimized microgrid and to be able to plot the trajectory

### study.jl

Set of functions to study the influence of artificial H2 price on the optimization by allowing to iterate over a range of values. 

    "function dataframe_result(result, H2_price_range, infinite_storage::Bool = false)"
and

    "function export_to_csv(df::DataFrame; prefix="result_for_H2_price")"

can be used to easily analyze the results with a CSV or directly with a dataframe.

## Changes TO DO

Implement multithreading when you iterate optimization for different H2 price as they are independent to one another.

File to modify :

    "study.jl"