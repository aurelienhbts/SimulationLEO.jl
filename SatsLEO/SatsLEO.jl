module SatsLEO

using LinearAlgebra, Plots

include("basis.jl")
export  μ, Re, ωe, Sat, deg2rad, rad2deg, R1, R3, 
        ecef_from_eci, latlon_from_ecef, eci_pos, walker_delta, myconstellation

include("calculations.jl")
export  visible, coverage_fraction, mean_coverage_fraction, eval_constellation

include("optimisation.jl")
export random_vec, mutate_vec, fitness, evolve_vec,
       mutate_vec_fixedN, fitness_fixedN, evolve_vec_fixedN,
       FITCACHE_GLOBAL, FITLOCK, FITCACHE_LOCAL, flush_fitcache!

include("visualisation.jl")
export show_coverage_heatmap, plot_constellation, plot_earth

end