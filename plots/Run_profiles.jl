"""
Plot the velocity profiles of the multigrid exercises and compare them with the analytical solution of the Shercliff flow.
"""

using DrWatson
@quickactivate "GridapMHDCalculations" 

using GridapMHDCalculations.plots: plot_Uz_profile
using GridapMHDCalculations: analytical_Hunt_u
using Plots
using Cubature

profile_dir_X=datadir("channel_mg_analysis", "ins_channel_L=4_b=1", "Ha=10_Re=1","vtk","profiles","Y_0")
profile_dir_X_2=datadir("channel_mumps_analysis", "ins_channel_L=4_b=1", "Ha=10_Re=1","vtk","profiles","Y_0")

filesnames_X = filter(endswith(".txt"), readdir(profile_dir_X))
filesnames_X_2 = filter(endswith(".txt"), readdir(profile_dir_X_2))

files_X = joinpath.(profile_dir_X, filesnames_X)
files_X_2 = joinpath.(profile_dir_X_2, filesnames_X_2)

push!(files_X, files_X_2...)

profile_dir_Y=datadir("channel_mg_analysis", "ins_channel_L=4_b=1", "Ha=10_Re=1","vtk","profiles","X_0")
profile_dir_Y_2=datadir("channel_mumps_analysis", "ins_channel_L=4_b=1", "Ha=10_Re=1","vtk","profiles","X_0")

filesnames_Y = filter(endswith(".txt"), readdir(profile_dir_Y))
filesnames_Y_2 = filter(endswith(".txt"), readdir(profile_dir_Y_2))

files_Y = joinpath.(profile_dir_Y, filesnames_Y)
files_Y_2 = joinpath.(profile_dir_Y_2, filesnames_Y_2)

push!(files_Y, files_Y_2...)

labels = [match(r"Nxy_c=.*?mg_levels=\d+", filename).match for filename in filesnames_X]
labels_2 = [match(r"Nxy=.*?mumps", filename).match for filename in filesnames_X_2]

push!(labels, labels_2...)

fig_side=plot_Uz_profile(files_X,"X";central_points_only=false,labels=String.(labels))

#I prefer to integrate directly compute the analycial flowrate becase at low Ha the kp of Shercliff is not that accurate 
Analytical_flowRate, error = hcubature(
    x -> analytical_Hunt_u(1.0,0.0,-1.0,10.0,100)(x[1],x[2]),
    [-1.0,-1.0],
    [1.0,1.0]
)

plot!(fig_side, 
      x->analytical_Hunt_u(1.0,0.0,-1.0,10.0,100)(x,0.0)*4/Analytical_flowRate,
      -1.0,1.0,
      label="Analytical",linewidth=1.5,
      linestyle=:dash
)
fig_Ha=plot_Uz_profile(files_Y,"Y",central_points_only=false,labels=String.(labels))

plot!(fig_Ha, 
      y->analytical_Hunt_u(1.0,0.0,-1.0,10.0,100)(0.0,y)*4/Analytical_flowRate,
      -1.0,1.0,
      label="Analytical",linewidth=1.5,
      linestyle=:dash
)

savefig(fig_side, "plots/figures/profiles/Uz_profile_X.png")
savefig(fig_Ha, "plots/figures/profiles/Uz_profile_Y.png")