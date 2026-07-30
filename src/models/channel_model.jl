"""
  channel_model(nc; kargs)
  
  )
  Function that returns an anonym function function that generates a Gridap model and the tags for BC. 
  Arguments are the geometrical characteristics and mesh characteristics of a rectangular cross-sectional channel.
#Arguments

-`nc: tuple of cells in each direction. 2D or 3D tuple selects a FD or 3D problem
-`b: channel aspect ratio
-`L: channel lenght
-`mesh_map = mesh map function
-`cw : wall conductance ratio (>100 --> perfect conductor)

"""
function channel_model(nc::Tuple{Int64,Int64,Int64};
    b = 1,
    L = 2,
    mesh_map = nothing,
    cw = 0
    )
    
    domain = (-b, b, -1.0, 1.0, 0.0, L)
    
    function (parts,rank_partition)  
    model=CartesianDiscreteModel(parts, rank_partition, domain, nc; map=mesh_map)
  
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
    	error("Thin wall BC not implemented in the SteadyState driver yet")
    end


    model, Dirichlet_Utags, Dirichlet_Jtags, Dirichlet_φtags
    end
end

function channel_model(nc::Tuple{Int64,Int64};
    b = 1,
    mesh_map = nothing
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

    model, Dirichlet_Utags, Dirichlet_Jtags
    end
end
