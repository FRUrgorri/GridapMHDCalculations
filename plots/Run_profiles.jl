using DrWatson
@quickactivate "GridapMHDCalculations" 

using GridapMHDCalculations.plots: plot_Uz_profile
using GridapMHDCalculations: analytical_Hunt_u, kp_shercliff_cartesian
using Plots

profile_dir_X=datadir("channel_mg_analysis", "ins_channel_L=4_b=1", "Ha=10_Re=1","vtk","profiles","Y_0")
filesnames_X = filter(endswith(".txt"), readdir(profile_dir_X))
files_X = joinpath.(profile_dir_X, filesnames_X)

profile_dir_Y=datadir("channel_mg_analysis", "ins_channel_L=4_b=1", "Ha=10_Re=1","vtk","profiles","X_0")
filesnames_Y = filter(endswith(".txt"), readdir(profile_dir_Y))
files_Y = joinpath.(profile_dir_Y, filesnames_Y)

fig_side=plot_Uz_profile(files_X,"X",central_points_only=true)

plot!(fig_side, 
      x->analytical_Hunt_u(1.0,0.0,-kp_shercliff_cartesian(1.0,10.0),10.0,100)(x,0.0),
      -1.0,1.0,
      label="Analytical",linewidth=1.5
)
fig_Ha=plot_Uz_profile(files_Y,"Y",central_points_only=true)

plot!(fig_Ha, 
      y->analytical_Hunt_u(1.0,0.0,-kp_shercliff_cartesian(1.0,10.0),10.0,100)(0.0,y),
      -1.0,1.0,
      label="Analytical",linewidth=1.5
)

savefig(fig_side, "plots/figures/profiles/Uz_profile_X.png")
savefig(fig_Ha, "plots/figures/profiles/Uz_profile_Y.png")