# Triangulation of the surface primitives, and the primitive list of a CSG tree.
#
# These are geometry, not plotting: a `Mesh` is what a primitive looks like as triangles, and both
# plotting backends consume it rather than producing it. They lived under `plotting/` alongside the
# Plots.jl recipes until those were removed (migration step 11), and moved here unchanged.

"""
    struct Mesh{T}

The triangulation of a surface primitive.

* `x`, `y`, `z`: vertex coordinates, in metres.
* `connections`: the vertices making up each face, as index lists.
"""
struct Mesh{T}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    connections::Vector{Vector{Int}}
end

"""
    mesh(p::AbstractSurfacePrimitive, n_arc)
    mesh(p::AbstractSurfacePrimitive, n_arc, n_vert_lines)

Triangulate `p`, dividing each full arc of a curved surface into `n_arc` segments.

The three-argument form additionally fixes how many lines run along the surface, which is what a
wireframe draws.
"""
function mesh(p::AbstractSurfacePrimitive{T}, n_arc::Int)::Mesh{T} where {T}
    vs = vertices(p, n_arc)
    x, y, z = broadcast(i -> getindex.(vs, i), (1,2,3))
    c = connections(p, n_arc)
    Mesh{T}(x,y,z,c)
end

function mesh(p::AbstractSurfacePrimitive{T}, n_arc::Int, n_vert_lines::Int)::Mesh{T} where {T}
    vs = vertices(p, n_arc)
    x, y, z = broadcast(i -> getindex.(vs, i), (1,2,3))
    c = connections(p, n_arc, n_vert_lines)
    Mesh{T}(x,y,z,c)
end

"""
    primitives(geometry) -> Tuple

The volume primitives a geometry is built from, flattened out of the CSG tree.
"""
@inline primitives(vp::AbstractVolumePrimitive) = (vp,)
@inline function primitives(csg::AbstractConstructiveGeometry)
    (
        primitives(csg.a)...,
        primitives(csg.b)...
    )
end

"""
    edges(e::Ellipse; n = 4) -> Vector{Edge}

Decompose an ellipse into `n` straight edges, which is how a curved outline is drawn.
"""
function edges(e::Ellipse{T,T,Nothing}; n = 4) where {T}
    φs = range(0, stop = 2π, length = n + 1)
    pts = [CartesianPoint(CylindricalPoint{T}(e.r, φ, zero(T))) for φ in φs]
    pts = map(p -> _transform_into_global_coordinate_system(p, e), pts)
    edges = [Edge(pts[i], pts[i+1]) for i in 1:n]
end
function edges(e::Ellipse{T,T,T}; n = 4) where {T}
    φs = range(T(0), stop = e.φ, length = n + 1)
    pts = [CartesianPoint(CylindricalPoint{T}(e.r, φ, zero(T))) for φ in φs]
    pts = map(p -> _transform_into_global_coordinate_system(p, e), pts)
    edges = [Edge(pts[i], pts[i+1]) for i in 1:n]
end
function edges(e::Ellipse{T,NTuple{2,Tuple{T}},Nothing}; n = 4) where {T}
    φs = range(0, stop = 2π, length = n + 1)
    pts = [CartesianPoint{T}(e.r[1][1] * cos(φ), e.r[2][1] * sin(φ), zero(T)) for φ in φs]
    pts = map(p -> _transform_into_global_coordinate_system(p, e), pts)
    edges = [Edge(pts[i], pts[i+1]) for i in 1:n]
end
