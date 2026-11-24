using Gmsh: gmsh

const R_Ha = 1.15
const R_side = 1.1
const N_Ha = 8
const N_side = 8
const n_inlet = 6
const n_outlet = 12
const r_exp = 1.0
const r = 4.0
const betta = 0.2
const L_in = 1.0
const L_out = 2.0

function add_physical_group(dim, name, entities)
    tag = gmsh.model.addPhysicalGroup(dim, entities)
    gmsh.model.setPhysicalName(dim, tag, name)
end

function set_progression_curve(signed_tag, num_nodes, ratio)
    coefficient = signed_tag < 0 ? 1 / ratio : ratio
    gmsh.model.mesh.setTransfiniteCurve(
        abs(signed_tag), num_nodes, "Progression", coefficient
    )
end

gmsh.initialize()
try
    gmsh.model.add("expansion")

    points = [
        (-L_in, 0.0),
        (-L_in, 1 / r),
        (0.0, 1 / r),
        (0.0, (1 + 1 / r) * 0.5),
        (0.0, 1.0),
        (L_out, 1.0),
        (L_out, (1 + 1 / r) * 0.5),
        (L_out, 1 / r),
        (L_out, 0.0),
        (L_out, -1 / r),
        (L_out, -(1 + 1 / r) * 0.5),
        (L_out, -1.0),
        (0.0, -1.0),
        (0.0, -(1 + 1 / r) * 0.5),
        (0.0, -1 / r),
        (-L_in, -1 / r),
        (0.0, 0.0),
    ]

    for (tag, (x, y)) in enumerate(points)
        gmsh.model.geo.addPoint(x, y, -betta, 1.0, tag)
    end

    lines = [
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7),
        (7, 8), (8, 9), (9, 10), (10, 11), (11, 12), (12, 13),
        (13, 14), (14, 15), (15, 16), (16, 1), (3, 17), (17, 15),
    ]

    for (tag, (start, finish)) in enumerate(lines)
        gmsh.model.geo.addLine(start, finish, tag)
    end

    gmsh.model.geo.addCurveLoop([1, 2, 17, 18, 15, 16], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addCurveLoop(
        [-17, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, -18], 2
    )
    gmsh.model.geo.addPlaneSurface([2], 2)

    # Build the extrusion before synchronizing the geometry to the model.
    N = N_side ÷ 2
    first_layer = 0.5 * (1 - R_side) / (1 - R_side^N)
    heights = Float64[]
    height = 0.0

    for i in 0:(N - 1)
        height += first_layer * R_side^i
        push!(heights, height)
    end
    for i in 0:(N - 1)
        height += first_layer * R_side^(N - i - 1)
        push!(heights, height)
    end

    gmsh.model.geo.extrude(
        [(2, 1), (2, 2)],
        0.0, 0.0, 2 * betta,
        ones(Int, N_side), heights, true,
    )

    gmsh.model.geo.synchronize()

    # Apply the mesh constraints after the final geometry synchronization.
    for curve in (-2, 15)
        set_progression_curve(curve, n_inlet, r_exp)
    end
    for curve in (5, -12)
        set_progression_curve(curve, n_outlet, r_exp)
    end
    for curve in (-1, 17, 8, 16, -18, -9)
        set_progression_curve(curve, N_Ha ÷ 2, R_Ha)
    end
    for curve in (3, -7, 13, -11, -14, 10, -4, 6)
        set_progression_curve(curve, N_Ha ÷ 2, R_Ha)
    end

    gmsh.model.mesh.setTransfiniteSurface(1, "Left", [2, 3, 15, 16])
    gmsh.model.mesh.setTransfiniteSurface(2, "Left", [5, 6, 12, 13])
    gmsh.model.mesh.setRecombine(2, 1)
    gmsh.model.mesh.setRecombine(2, 2)

    add_physical_group(2, "inlet", [29, 49])
    add_physical_group(2, "outlet", [85, 89, 93, 97, 101, 105])
    add_physical_group(2, "wall", [1, 2, 33, 45, 109, 81, 50, 122, 77, 73, 117, 113])
    add_physical_group(3, "PbLi", [1, 2, 3, 4])

    add_physical_group(1, "inlet", [27])
    add_physical_group(1, "outlet", [84, 88, 92, 96, 100])
    add_physical_group(1, "wall", [
        1, 16, 20, 25, 28, 44, 2, 15, 21, 24, 3, 4, 17, 18, 13, 14,
        53, 54, 22, 23, 63, 64, 76, 72, 32, 40, 112, 108, 5, 12, 55,
        62, 6, 7, 8, 9, 10, 11, 56, 57, 58, 59, 60, 61, 80, 104,
    ])
    
    add_physical_group(0, "wall", [
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 23, 27, 31, 35, 41, 45, 49, 53, 57, 61, 65, 69, 73,
        77, 81,
    ])

    gmsh.model.mesh.generate(3)
    gmsh.write(joinpath(@__DIR__, "expansion.msh"))
finally
    gmsh.finalize()
end