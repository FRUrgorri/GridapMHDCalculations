#Velocity Profiles
function u_parabolic(b::Real)
    (x,y) -> (9/(4*b^2))*(x^2 - b^2)*(y^2 - 1)
end

"""
  analytical_GeneralHunt_u(l, d_b, grad_p, Ha, n, x)

General Hunt analytical formula (d_b = 0 for Shercliff).

# Arguments
- `l::Float64`: channel aspect ratio.
- `d_b::Float64`: Hartmann walls conductivity ratio.
- `grad_p::Float64`: dimensionless (MHD version) presure gradient.
- `Ha::Float64`: Hartmann number.
- `n::Int`: number of sumands included in Fourier series.
"""
function analytical_Hunt_u(
  # General Hunt analytical formula (d_b = 0 for Shercliff)
  l::Real,          # channel aspect ratio
  d_b::Real,        # Hartmann walls conductivity ratio
  grad_p::Real,     # Dimensionless (MHD version) presure gradient
  Ha::Real,         # Hartmann number
  n::Int           # number of sumands included in Fourier series
 )
  
  function V(x::Real,y::Real)
    _V = 0.0
    for k in 0:n
      α_k = (k + 0.5)*π/l
      N = (Ha^2 + 4*α_k^2)^(0.5)
      r1_k = 0.5*( Ha + N) #Typo in the Ni. et al. (2013) paper, it is N not N^2
      r2_k = 0.5*(-Ha + N)

      eplus_1k = 1 + exp(-2*r1_k)
      eminus_1k = 1 - exp(-2*r1_k)
      eplus_2k = 1 + exp(-2*r2_k)
      eminus_2k = 1 - exp(-2*r2_k)
      eplus_k = 1 + exp(-2*(r1_k+r2_k))
    
      e_x_1k = 0.5*(exp(-r1_k*(1-y))+exp(-r1_k*(1+y)))
      e_x_2k = 0.5*(exp(-r2_k*(1-y))+exp(-r2_k*(1+y)))

      V2 = ((d_b*r2_k + eminus_2k/eplus_2k)*e_x_1k)/(0.5*N*d_b*eplus_1k + eplus_k/eplus_2k)
      V3 = ((d_b*r1_k + eminus_1k/eplus_1k)*e_x_2k)/(0.5*N*d_b*eplus_2k + eplus_k/eplus_1k)

      _V += 2*(-1)^k*cos(α_k * x)/(l*α_k^3) * (1-V2-V3)
    end
    return _V
  end

  return (x,y) -> V(x,y)*Ha^2*(-grad_p)
end

#Analitical formulas for pressure drop gradients (:mhd discretization of the pressure field)
kp_shercliff_cartesian(b::Real,Ha::Real)= 1/(Ha*(1-0.852*Ha^(-0.5)/b-1/Ha))
kp_shercliff_cylinder(Ha::Real) = (3/8)*pi/(Ha-(3/2)*pi)
kp_hunt(b::Real,Ha::Real) = 1/(Ha*(1-0.956*Ha^(-0.5)/b-1/Ha))
kp_glukhih(Ha::Real,cw::Real) = (3/8)*pi*(1+0.833*cw*Ha-0.019*(cw*Ha)^2)/Ha

function kp_tillac(b::Real,Ha::Real,cw_s::Real,cw_Ha::Real)
  k_s = (1/(3*b))*(Ha^(0.5)/(1+cw_s*Ha^(0.5)))
  k_Ha = (1+cw_Ha)/(1/Ha + cw_Ha)
  return 1/(k_s+k_Ha)
end
 