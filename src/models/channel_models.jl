"""
Abstract supertype for all the channel models to be run with this package
"""
abstract type channel_models <:mounted_models end

"""
 Struct containing the channel geometry fields (in a way equivalent to Gridap CartesianDescriptor but simpler)

 #Fields
    `b:` Channel aspect ratio (x direction)
    `L:` Channel lenght (z direction)
"""

struct channel_geom
    b::Real
    L::Real
end

(ch::channel_geom)() = (-ch.b, ch.b, -1.0, 1.0, 0.0, ch.L)

"""
 Struct containing the mesh information for a channel

 #Fields
    `nc: Tuple with the number of cells of the coarser multigrid level (or singlegrid level)`
    `map:` map function to be applied to the mesh
    `level:` Number of multigrid levels
    `nrefs:` Refinement factor per multigrid level
"""

struct channel_mesh
    nc::NTuple{3,Integer}
    map::Function
    levels::Integer;
    nrefs::Union{Integer,NTuple{3,Integer}} 
end

channel_mesh(nc,map) = channel_mesh(nc,map,1,1)
channel_mesh(nc) = channel_mesh(nc,identity)


"""
Struct containign the information of a GridapMHD insulated channel model
#Fields
 -`b`: channel aspect ratio
 -`L`: channel length
 -`nc`: Number of cells accross each direction of the channel 
 -`BCs`: tags and BC values of the channel model
 -`B`: External magnetic field
 -`source`: external source
"""

struct insulated_channel <: channel_models
    geom::channel_geom
    mesh::channel_mesh
    BCs::BC
    B::Union{VectorValue{3,Float64},Function}          
    source::Union{VectorValue{3,Float64},Function}
end

insulated_channel(geom,mesh,BCs,B) = insulated_channel(geom,mesh,BCs,B,VectorValue(0.0,0.0,0.0))
insulated_channel(geom,mesh,BCs) = insulated_channel(geom,mesh,BCs,VectorValue(0.0,1.0,0.0))

"""
Instance of insulated channel model in serial
"""
function (ins_ch::insulated_channel)()
    domain = ins_ch.geom()
    nc = ins_ch.mesh.nc
    map = ins_ch.mesh.map
    model=CartesianDiscreteModel(domain, nc, map)
    add_insulated_tags!(model, ins_ch.BCs.tags)
    return model
end

"""
Instance of insulated channel model with mpi or sequential(for debugging) backends
"""
function (ins_ch::insulated_channel)(parts::Union{AbstractVector,MPIArray},ranks::NTuple{3,Integer})
    
    domain = ins_ch.geom()
    nc = ins_ch.mesh.nc
    map = ins_ch.mesh.map

    if ins_ch.mesh.levels == 1 #Single grid
       
        model=CartesianDiscreteModel(parts, ranks, domain, nc,map)
        add_insulated_tags!(model, ins_ch.BCs.tags)
        multigrid = nothing

    else #multigrid

        ranks_per_level = fill(ranks,ins_ch.mesh.levels) #Same amount of processors per level (potentially troublesome in coarse levels)
        model_hierarchy = CartesianModelHierarchy(parts, ranks_per_level, domain, nc; map=map, nrefs = ins_ch.mesh.nrefs)
        multigrid = Dict{Symbol,Any}(
         :mh => model_hierarchy,
         :num_refs_coarse => 0,  #What is this?
         :ranks_per_level => ranks_per_level,
        ) 
        add_insulated_tags!(model_hierarchy, ins_ch.BCs.tags) 
        model = get_model(model_hierarchy,1)
        
    end
    
    return model, multigrid
end

"""
Function for adding the insulated tags to a Gridap channel model
    
#Arguments

-`model: Gridap model 
-`tags::BC_tags: tags to be added 
"""

function add_insulated_tags!(model::Union{CartesianDiscreteModel,GridapDistributed.DistributedDiscreteModel},tags::BC_tags)
    labels = get_face_labeling(model)
    tags_inlet = append!(collect(1:4), [9, 10, 13, 14], [21])
    tags_outlet = append!(collect(5:8), [11, 12, 15, 16], [22])
    tags_wall = append!(collect(1:20), [23, 24, 25, 26])

    #Use the names in J tag since everything is insulated
    add_tag_from_tags!(labels, tags.tags_J[1], tags_inlet)
    add_tag_from_tags!(labels, tags.tags_J[2], tags_outlet)
    add_tag_from_tags!(labels, tags.tags_J[3], tags_wall)

    return nothing
end

#Add tag to every level
function add_insulated_tags!(mh::ModelHierarchy,tags::BC_tags)
    map(mh) do mh_level
        m_level = get_model(mh_level)
        add_insulated_tags!(m_level,tags)
    end
    return nothing
end