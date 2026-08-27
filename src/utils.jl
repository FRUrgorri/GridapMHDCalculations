function u_parabolic(b::Real)
    u(x,y) = (9/(4*b^2))*(x^2 - b^2)*(y^2 - 1)
end