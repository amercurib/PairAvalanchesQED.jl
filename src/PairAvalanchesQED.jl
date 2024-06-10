module PairAvalanchesQED
#export path_table_b0
#path_table_b0 = joinpath(pkgdir(PairAvalanchesQED),"tables","b0_table.jld")
#export path_table_c0 
#path_table_c0= joinpath(pkgdir(PairAvalanchesQED),"tables","b0_table.jld")

#SFQED Rates part
export b0_low_chi, b0_high_chi, b0, Wnbw_SI, Wnbw, b0_from_integral, compute_b0_table, c0_low_chi, c0_high_chi, c0, Wncs_SI, Wncs, c0_from_integral,  compute_c0_table
include("Rates.jl")
using .Rates

#Short time dynamics part
export chi_e_SI, gamma_e_SI, chi_e, gamma_e, Wgam_t_SI, Wpair_t_SI, Wgam_t, Wpair_t
include("ShortDynamics.jl")
using .ShortDynamics


#Avalanche model part
export int_Wncs, solve_tem, chi_em, gamma_em, GR_steady_state, GR_constant_MidIntensity, GR_pure_rotating_MidIntensity, GR_constant_MidIntensity_eps2weff, GR_constant_HighIntensity, GR_pure_rotating_HighIntensity, GR_constant_HighIntensity_eps2weff, GR_constant, GR_pure_rotating, GR_constant_eps2weff
include("AvalancheModelRotatingConfigs.jl")
using .AvalancheModelRotatingConfigs

export computing_invariants, make_TF, make_square_field_tensor, computing_derivatives_EB, make_derivative_tensor, make_J_matrix, invert_J_matrix, computing_eigenvector_f1, lower_f1_index, compute_weff, weff_epsilon_from_field
#Computation of invariant and quantities for short time dynamics from a field configuration
include("EpsilonOmegaEff.jl")
using .EpsilonOmegaEff

end # module PairAvalanchesQED
