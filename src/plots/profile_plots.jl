"""
Functions for plotting profiles from CSV or txt files generated with paraview
    #Arguments
        files::Vector{String}: Vector of paths to the profile files
        first_header::String: Name of the first column in the profile file
        second_header::String: Name of the second column in the profile file
"""
function read_profile_columns(file::String, first_header, second_header)
    df = CSV.read(file, DataFrame)

    first_col = Symbol(first_header)
    second_col = Symbol(second_header)

    first_col in propertynames(df) ||
        throw(ArgumentError("Column \"$first_header\" was not found in $file"))
    second_col in propertynames(df) ||
        throw(ArgumentError("Column \"$second_header\" was not found in $file"))

    return df[!, first_col], df[!, second_col]
end

function read_profile_columns(files::Vector{String}, first_header, second_header)
    return map(file -> read_profile_columns(file, first_header, second_header), files)
end

function _central_step_indices(position, uz)
    length(position) == length(uz) ||
        throw(DimensionMismatch("Position and Uz columns must have the same length"))

    central_indices = Int[]
    run_start = 1
    while run_start <= length(uz)
        run_end = run_start
        while run_end < length(uz) && isequal(uz[run_end + 1], uz[run_start])
            run_end += 1
        end

        step_positions = view(position, run_start:run_end)
        midpoint = (minimum(step_positions) + maximum(step_positions)) / 2
        central_index = run_start - 1 + argmin(abs.(step_positions .- midpoint))
        push!(central_indices, central_index)
        run_start = run_end + 1
    end

    return central_indices
end


"""
Function to plot Uz profiles from CSV or txt files generated with paraview
    #Arguments
        files::Vector{String}: Vector of paths to the profile files
        side::String: Side of the channel where the profile was taken, either "X" or "Y"
    #Keyword arguments
        central_points_only::Bool: Plot only the point nearest the center of each constant-Uz step
"""
function plot_Uz_profile(
    files::Vector{String},
    side::String;
    central_points_only::Bool=false,
    labels::Vector{String}=["" for _ in 1:length(files)],
)
    isempty(files) && throw(ArgumentError("At least one profile file is required"))

    side ∈ ("X", "Y") || throw(ArgumentError("Side must be either \"X\" or \"Y\""))

    side == "X" ? xlabel="X" : xlabel ="Y"
    side == "X" ? header = "Points:0" : header = "Points:1"

    p = plot(xlabel=xlabel, ylabel="Uz", legend=:bottom)

    for ((i,), (position, uz)) in zip(
        enumerate(files),
        read_profile_columns(files, header, "uh:2"),
    )

        if central_points_only
            central_indices = _central_step_indices(position, uz)
            plot!(
                p,
                position[central_indices],
                uz[central_indices];
                label=labels[i],
                linewidth=1.5
            )
        else
            plot!(p, position, uz; label=labels[i], linewidth=1.5)
        end
    end

    return p
end
