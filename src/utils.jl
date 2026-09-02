#Velocity Profiles
function u_parabolic(b::Real)
    u(x,y) = (9/(4*b^2))*(x^2 - b^2)*(y^2 - 1)
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
 