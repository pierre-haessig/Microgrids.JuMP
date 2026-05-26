# Microgrids.JuMP

Code for microgrid sizing optimization using an algebraic formulation ([implemented with JuMP](https://jump.dev/)), rather than a simulator-based formulation (e.g. with [Microgrids.jl](https://github.com/Microgrids-X/Microgrids.jl)).

It is mostly a classical/simplified Linear Programming based energy system optimization model like many others (DER-CAM, OSeMOSYS, REopt, PyPSA, FINE, MicroGridsPy, MANGO, Tulipa). Still, it comes with two nice "twists":

- Available as a self-contained literate programming version as a **single Jupyter notebook** to make it easy to grasp for fresh researchers → [Microgrid_optimization_JuMP.ipynb](Microgrid_optimization_JuMP.ipynb)
  - (but also available as a set of encapsulated functions in [Microgrid_JuMP_common.jl](Microgrid_JuMP_common.jl) for reuse in numerical experiments, like sensitivity analysis, model comparison...)
- Implements a **convex system cost model** which allows **usage-dependent component lifetimes**, along its piecewise linear implementation (see explanations in the manuscript referenced below)

## Installation and dependencies

Three Julia packages need to be installed to use Microgrids.JuMP

- [Microgrids.jl](https://github.com/Microgrids-X/Microgrids.jl): to hold microgrid parameters such as prices, load curve...
- [JuMP](https://jump.dev/): for the optimization model formulation
- a Linear Programming solver, for example [HiGHS.jl](github.com/jump-dev/HiGHS.jl)

Also, most code is available as  Jupyter notebooks, so [IJulia.jl](https://github.com/JuliaLang/IJulia.jl) is needed as well to interact with those.

As of now, Microgrids.JuMP is not a Julia package, so it's installation instructions are: 

1. clone/copy this repo
1. use the code, starting with the notebook [Microgrid_optimization_JuMP.ipynb](Microgrid_optimization_JuMP.ipynb)

## Academic use and theory

### Article on convex cost with usage-dependent component lifetimes

This code is the implementation of the numerical results of the yet-to-be-published manuscript *“Convexity of the annualized cost of energy systems with usage-dependent component lifetimes”* https://hal.science/hal-05052473.

This article gives the lifetime model description, convexity proof and a generic piecewise linear implementation. 

The results of this article are generated in the notebook [Cvx_article_experiments.ipynb](Cvx_article_experiments.ipynb)

### Talk on the anticipativity of optimization-based energy managment with perfect foresight

The first version of this code was used in a [flash talk about anticipativity](https://forum.openmod.org/t/grenoble-workshop-2024-lightning-talks/4438/34) at the [Openmod workshop 2024 in Grenoble](https://forum.openmod.org/t/grenoble-workshop-2024-registration-and-details/4350).

The results presented in the talk are generated in the notebook [Microgrid_JuMP-BB_comparison.ipynb](Microgrid_JuMP-BB_comparison.ipynb)

### Infinite storage model

The [infinite_storage](https://github.com/pierre-haessig/Microgrids.JuMP/tree/infinite_storage) branch, developed with https://github.com/Goussedevanilles, contains an experimental model variant where a finite size hydrogen storage is replaced by an infinite storage (never full nor empty) where the usage of H2 is penalized at a fixed price.
