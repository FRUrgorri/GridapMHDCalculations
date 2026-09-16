"""
  channel_model(nc::NTuple{3,<:Integer}; kargs...)
  
  Function that returns an anonym function function that generates a Gridap model and the tags for BC. 
  Arguments are the geometrical characteristics and mesh characteristics of a rectangular cross-sectional channel.

  #Arguments

    -`nc: tuple of cells in each direction. 2D or 3D tuple selects a FD or 3D problem
    -`b: channel aspect ratio
    -`L: channel lenght
    -`mesh_map = mesh map function
    -`cw : wall conductance ratio (>100 --> perfect conductor)

"""
function channel_model(nc::NTuple{3,Integer};
    b::Real = 1.0,
    L::Real = 2.0,
    mesh_map::Function = identity,
    cw::Real = 0.0
    )
    
    domain = (-b, b, -1.0, 1.0, 0.0, L)
    
    function (parts::AbstractVector{Integer},rank_partition::Union{Integer,NTuple{3,Integer}})  
        
        model=CartesianDiscreteModel(parts, rank_partition, domain, nc; map=mesh_map)
  
        tags = add_channel_tags!(model;cw=cw)

        model, tags, Nothing
    end
end

#Method for a Cartesian Hierarchy model
function channel_model(nc::NTuple{3,Integer}, levels::Integer;
    nrefs::Union{<:Integer,NTuple{3,Integer}} = 2,
    b::Real = 1.0,
    L::Real = 2.0,
    mesh_map::Function = identity,
    cw::Real = 0.0,
    )
    
    domain = (-b, b, -1.0, 1.0, 0.0, L)

    function (parts::AbstractVector{Integer},rank_partition::NTuple{3,Integer})  
    
        ranks_per_level = fill(rank_partition,levels) #Same amount of processors per level (potentially troublesome in coarse levels)

        m_hierarchy = CartesianModelHierarchy(parts, ranks_per_level, domain, nc; map=mesh_map, nrefs=nrefs)      #This has no redistribution (see GridapDistributed -> Redistribution.jl)

        multigrid = Dict{Symbol,Any}(
         :mh => m_hierarchy,
         :num_refs_coarse => 0,  #What is this?
         :ranks_per_level => ranks_per_level,
        )

        tags = add_channel_tags!(m_hierarchy;cw=cw)

        model = get_model(m_hierarchy,1)     #Pass the finnest level of the hierarchy as model (equivalent to get_model(mh[1]))
    
        model, tags, multigrid
    end
end

"""
add_channel_tags!(model::Union{CartesianDiscreteModel,ModelHierarchyLevel},cw::Real)

Add the tags "inlet", "outlet" "insulating" and "conducting" to the tag list of the model face labeling
Returns arrays of the Dirichlet tags for U, J and Phi so they can be used for the BC definition

"""
function add_channel_tags!(model::Union{CartesianDiscreteModel,GridapDistributed.DistributedDiscreteModel};cw::Real = 0.0)
    # Vertex tags: [1:8]
    # Edge tags: [9:20]
    # Surf tags: [21:26]


    labels = get_face_labeling(model)
    tags_inlet = append!(collect(1:4), [9, 10, 13, 14], [21])
    tags_outlet = append!(collect(5:8), [11, 12, 15, 16], [22])
    tags_insulated = append!(collect(1:20), [23, 24, 25, 26])
    add_tag_from_tags!(labels, "inlet", tags_inlet)
    add_tag_from_tags!(labels, "outlet", tags_outlet)
    
    if iszero(cw)
    	add_tag_from_tags!(labels, "insulated", tags_insulated)
	
    	#Neumann tags are the default, so there is no need to specified "outlet" as Neumann for example 
    	Dirichlet_Utags=["inlet","insulated"]
    	Dirichlet_Jtags=["inlet","outlet","insulated"]
    	Dirichlet_φtags=[]
    elseif cw > 100 
    	add_tag_from_tags!(labels, "conducting", tags_insulated)
	    Dirichlet_Utags=["inlet","conducting"]
        Dirichlet_Jtags=["inlet","outlet"]
        Dirichlet_φtags=["conducting"]
    else
    	error("Thin wall BC not implemented in the SteadyState driver yet")  #TBD
    end

    return Dirichlet_tag_names(Dirichlet_Utags, Dirichlet_Jtags, Dirichlet_φtags)
end

#Collect tag every level
function add_channel_tags!(mh::ModelHierarchy;cw=cw)
    map(mh[1:end-1]) do mh_level
        m_level = get_model(mh_level)
        add_channel_tags!(m_level;cw=cw)
    end
        m_level_last = get_model(mh[end])
        tags = add_channel_tags!(m_level_last;cw=cw)
    return tags
end

"""
#2D method for FD calculations (To be updated)
function channel_model(nc::NTuple{2,<:Integer};
    b::Real = 1.0,
    mesh_map::Function = identity
    )
    
    domain = (-b, b, -1.0, 1.0, 0.0, 0.1)
    _nc = (nc[1],nc[2],3)

    function (parts,rank_partition)

    model = CartesianDiscreteModel(parts, rank_partition, domain, _nc; isperiodic=(false,false,true),  map=mesh_map)
   
    labels = get_face_labeling(model)
    tags_insulated = append!(collect(1:20), [23, 24, 25, 26])
    add_tag_from_tags!(labels, "insulated", tags_insulated)

    Dirichlet_Utags=["insulated"]
    Dirichlet_Jtags=["insulated"]

    tags = (Dirichlet_Utags, Dirichlet_Jtags)

    model, tags
    end
end
"""
