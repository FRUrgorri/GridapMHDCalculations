function collect_and_analyze_results(analysis_dir::String)

    monitor_results = collect_results(
        joinpath(analysis_dir, "monitors");
        valid_filetypes=".bson",
    )
    info_results = collect_results(
        joinpath(analysis_dir, "vtk", "info");
        valid_filetypes=".bson",
    )

    #Sort the results by path and check that they contain the same filenames. Then, merge the two DataFrames into a single DataFrame.
    sort!(monitor_results, :path)
    sort!(info_results, :path)
    basename.(monitor_results.path) == basename.(info_results.path) ||
        throw(ArgumentError("The two result sets do not contain the same filenames"))

    duplicate_names = intersect(names(monitor_results), names(info_results))
    select!(info_results, Not(duplicate_names))
    results = hcat(monitor_results, info_results)

    results[!, :Nxy_int] = results.Nxy_c .* 2 .^ (results.mg_levels .- 1) #Create a new column that is the product of Nxy_c and mg_levels (the integration mesh)
    results[!, :mapFunction_name] = [                                     #Create a new column that contains the name of the mapFunction for each row.
        map_function == identity ? "identity" : "Roberts"
        for map_function in results.mapFunction
    ]

    sort!(results, :Nxy_int)        #Create a dataframe group for each unique value of Nxy_int
    results_by_Nxy_int = groupby(results, :Nxy_int) 

    #Filter the results to only include those that converged (ls1_residuals < 0.001)
    converged_mask = last.(results.ls1_residuals) .< 0.001
    converged = results[converged_mask, :]
    not_converged = results[.!converged_mask, :]

    converged[!, :Nxy_int] = converged.Nxy_c .* 2 .^ (converged.mg_levels .- 1) #Create a new column that is the product of Nxy_c and mg_levels (the integration mesh)
    converged[!, :mapFunction_name] = [                                     #Create a new column that contains the name of the mapFunction for each row.
       map_function == identity ? "identity" : "Roberts"
       for map_function in converged.mapFunction
    ]

    sort!(converged, :Nxy_int)        #Create a dataframe group for each unique value of Nxy_int
    converged_by_Nxy_int = groupby(converged, :Nxy_int) 

    return results_by_Nxy_int, converged_by_Nxy_int, not_converged
end

function plot_Uz_vs_μ(
    df_g::GroupedDataFrame;
    map_function_name::Union{Nothing,String}=nothing,
)
    map_function_name in (nothing, "identity", "Roberts") ||
        throw(ArgumentError("map_function_name must be nothing, \"identity\", or \"Roberts\""))

    nplots = length(df_g)
    ncols = ceil(Int, sqrt(nplots))
    nrows = ceil(Int, nplots / ncols)
    uz_wall_plot = plot(
        layout=(nrows, ncols),
        size=(500 * ncols, 350 * nrows),
    )
    map_colors = Dict("identity" => :blue, "Roberts" => :orange)
    line_styles = [:solid, :dash, :dot, :dashdot]
    marker_shapes = [:circle, :rect, :diamond, :utriangle]

    for (plot_index, mesh_results) in enumerate(df_g)
        plot!(
            uz_wall_plot[plot_index];
            title="Nxy_int = $(first(mesh_results.Nxy_int))",
            xlabel="μ",
            ylabel="Uz wall",
        )

        for map_results in groupby(mesh_results, :mapFunction_name)
            map_name = first(map_results.mapFunction_name)
            map_function_name === nothing || map_name == map_function_name || continue

            level_groups = sort(
                collect(groupby(map_results, :mg_levels));
                by=group -> first(group.mg_levels),
            )
            for (level_index, level_results) in enumerate(level_groups)
                curve = sort(level_results, :μ_BC)
                plot!(
                    uz_wall_plot[plot_index],
                    curve.μ_BC,
                    curve.Uz_wall;
                    label="$map_name, mg_levels=$(first(curve.mg_levels))",
                    color=get(map_colors, map_name, :gray),
                    linestyle=line_styles[mod1(level_index, length(line_styles))],
                    marker=marker_shapes[mod1(level_index, length(marker_shapes))],
                )
            end
        end
    end

    return uz_wall_plot
end
