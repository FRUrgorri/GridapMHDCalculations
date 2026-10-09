using DrWatson
@quickactivate "GridapMHDCalculations"

using GridapMHDCalculations.plots: collect_and_analyze_results, plot_Uz_vs_μ
using Plots

analysis_dir = datadir("channel_mg_analysis", "ins_channel_L=4_b=1", "Ha=10_Re=1")
results_by_Nxy_int, converged_by_Nxy_int, not_converged = collect_and_analyze_results(analysis_dir)

Full_Uz_wall_plot = plot_Uz_vs_μ(results_by_Nxy_int)
Converged_Uz_wall_plot = plot_Uz_vs_μ(converged_by_Nxy_int)
Converged_Uz_wall_Roberts_plot = plot_Uz_vs_μ(converged_by_Nxy_int; map_function_name="Roberts")

plotsdir =projectdir("plots/figures")
savefig(Full_Uz_wall_plot, joinpath(plotsdir, "Full_Uz_wall_plot.png"))
savefig(Converged_Uz_wall_plot, joinpath(plotsdir, "Converged_Uz_wall_plot.png"))
savefig(Converged_Uz_wall_Roberts_plot, joinpath(plotsdir, "Converged_Uz_wall_Roberts_plot.png"))