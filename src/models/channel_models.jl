"""
Abstract supertype for all the channel models to be run with this package
"""
abstract type channel_models <:mounted_models end

struct channel_geom
    b::Real
    L::Real
end

(ch::channel_geom)() = (-ch.b, ch.b, -1.0, 1.0, 0.0, ch.L)

struct channel_mesh
    nc::NTuple{3,Integer}
    map::Function
end

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
    model=CartesianDiscreteModel(parts, ranks, domain, nc,map)
    add_insulated_tags!(model, ins_ch.BCs.tags)
    return model
end

"""
struct insulated_mg_channel <: channel_models
    b::Real
    L::Real
    nc::NTuple{3,Integer}
    levels::Integer
    BCs::BC
    B::Union{VectorValue{3,Float64},Function}
    source::Union{VectorValue{3,Float64},Function}
end

insulated_channel(b,L,nc,levels,BCs,B) = insulated_channel(b,L,nc,levels,BCs,B,VectorValue(0.0,0.0,0.0))
"""

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

