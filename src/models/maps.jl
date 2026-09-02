"""
  Mesh streching according to the Roberts formula described in S. Smolentsev et al. (2005) 10.1016/j.fusengdes.2005.01.003
"""
function stretchMHD(coord;domain=(0.0,1.0,0.0,1.0,0.0,1.0),factor=(1.0,1.0,1.0),dirs=(1,2,3))
  ncoord = collect(coord.data)
  for (i,dir) in enumerate(dirs)
    ξ0 = domain[i*2-1]
    ξ1 = domain[i*2]
    l =  ξ1 - ξ0
    c = (factor[i] + 1)/(factor[i] - 1)

    if l > 0
      if ξ0 <= coord[dir] <= ξ1
        ξx = (coord[dir] - ξ0)/l                     # ξx from 0 to 1 uniformly distributed
        ξx_streched = factor[i]*(c^ξx-1)/(1+c^ξx)    # ξx streched from 0 to 1 towards 1
        ncoord[dir] =  ξx_streched*l + ξ0            # coords streched towards ξ1
      end
    else
      if ξ1 <= coord[dir] <= ξ0
        ξx = (coord[dir] - ξ0)/l                     # ξx from 0 to 1 uniformly distributed
        ξx_streched = factor[i]*(c^ξx-1)/(1+c^ξx)    # ξx streched from 0 to 1 towards 1
        ncoord[dir] =  ξx_streched*l + ξ0            # coords streched towards ξ1
      end
    end
  end
  return VectorValue(ncoord)
end


"""
    map_Roberts(b,Ha)
    
Function that returns a function for the cross sectional map using Roberts formula.
See strechMHD function for detailed formula

# Arguments
- `b` :  Channel aspect ratio
- `Ha: Streching factor along the direction perpendicular to B
    
"""
function map_Roberts(b,Ha)
     
     stretch_Ha = sqrt(Ha/(Ha-1))
     stretch_side = sqrt(sqrt(Ha)/(sqrt(Ha)-1))
     
     function (coord)
       ncoord = stretchMHD(coord,domain=(0,-b,0,-1.0),factor=(stretch_side,stretch_Ha),dirs=(1,2))
       ncoord = stretchMHD(ncoord,domain=(0,b,0,1.0),factor=(stretch_side,stretch_Ha),dirs=(1,2))
       ncoord
     end
end

