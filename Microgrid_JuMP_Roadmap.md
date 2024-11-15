# ROADMAP

## Code change TODOs, identified on 15 Nov 2024:

- move all parameters in `model_data` (already WIP)
- use `mg_base` to store the maximat ratings (instead of global variables)
- remove the need for global function `create_microgrid`  in `setup_optim_mg_jump` from notebooks (or internalize to `Microgrid_JuMP_common.jl`)



## Code outline 



*version sept 2024* → outdated



### JuMP optim

```
function build_optim_mg_stage!(mg, model_data::Dict{String,Any};
        shed_max = 0.0,
        ndays=365,
        fixed_lifetimes=false,
        relax_gain = 1.0,
    )
```

1. create_microgrid 2x
1. create JuMP model 2x
1. calls build_optim_mg_stage! 2x
1. sets objective

```    
function setup_optim_mg_jump(optimizer;
        shed_max=0.0,
        ndays=365,
        fixed_lifetimes=false,
        relax_gain = 1.0,
        create_mg_base=create_mg_base)
```

1. create_microgrid
2. create JuMP model
3. calls `build_optim_mg_stage!`
4. sets objective
```
function optim_mg_jump(optimizer;
        shed_max=0.0,
        ndays=365,
        fixed_lifetimes=false,
        relax_gain = 1.0,
        create_mg_base=create_mg_base,
        model_custom=nothing  # used to force sizing
    )
```
1. calls `setup_optim_mg_jump`

2. calls model_custom is any
3. solve
4. extract results + diagnostics

### Comparison

function optim_mg_both(shed_max;
        jump_optimizer,
        relax_gain = 1.0,
        bb_algo=:GN_CRS2_LM, bb_srand=1,
        create_mg_base = create_mg_modified,
    )
    
    calls optim_mg_jump
    calls optim_mg_bb (and then obj_multi)
    
    and then simulate x_alg

### Two stage

```julia
function setup_optim_mg_jump_2SP(shed_max, optimizer, Tmid;
        ndays=365,
        fixed_lifetimes=false,
        relax_gain = 1.0,
        investment_price_sto=investment_price_sto,
        investment_price_gen=investment_price_gen,
        fuel_price1=fuel_price,
        fuel_price_scenar=fuel_price_scenar,
        proba_scenar=proba_scenar
    )
```

```julia
function optim_mg_jump_2SP(shed_max, optimizer, Tmid;
        model_custom=nothing,
        ndays=365,
        fixed_lifetimes=false,
        relax_gain=1.0,
        investment_price_sto=investment_price_sto,
        investment_price_gen=investment_price_gen,
        fuel_price1=fuel_price,
        fuel_price_scenar=fuel_price_scenar,
        proba_scenar=proba_scenar
    )
```

