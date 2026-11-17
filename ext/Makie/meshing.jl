# SSD meshes to Makie geometry.
#
# `ConstructiveSolidGeometry.Mesh{T}` carries three coordinate vectors and a `connections` vector
# whose entries mean two different things depending on which `mesh` method produced it (plan §1.8):
#
#   * `mesh(s, n_arc)` gives one entry per **face**, listing the vertex indices around it. Faces are
#     quads on the mantles, triangles on a circular area and N-gons on the flat end of a prism, so
#     `QuadFace` alone is not enough and every face is fanned into triangles here.
#   * `mesh(s, n_arc, n_vert_lines)` gives one entry per **line**, a pair of vertex indices.
#
# The Plots recipes fed both to the same `:mesh3d` seriestype and let GR sort it out. Makie has the
# two idiomatic calls instead: `mesh!` on a `GeometryBasics.Mesh` for the first, `linesegments!` on a
# flat point vector for the second.
#
# Coordinates are rescaled to display units on the way out (D7), like every other recipe here.

"""
    mesh_points(m::CSG.Mesh) -> Vector{Point3f}

Give the vertices of `m` in display units.
"""
function mesh_points(m::Mesh)
    u = internal_length_unit
    x, y, z = display_values(m.x, u), display_values(m.y, u), display_values(m.z, u)
    Makie.Point3f[Makie.Point3f(x[i], y[i], z[i]) for i in eachindex(x)]
end

"""
    triangle_faces(connections) -> Vector{GLTriangleFace}

Triangulate face-index lists by fanning each face from its first vertex.

Valid because every face SSD's meshing produces is planar and convex; a fan across a non-convex face
would cover area outside it.
"""
function triangle_faces(connections)
    faces = Makie.GeometryBasics.GLTriangleFace[]
    for f in connections
        length(f) < 3 && continue
        for i in 2:(length(f) - 1)
            push!(faces, Makie.GeometryBasics.GLTriangleFace(f[1], f[i], f[i + 1]))
        end
    end
    faces
end

"""
    MeshBuffer()

Accumulator that concatenates several `CSG.Mesh`es into one `GeometryBasics.Mesh`.

One `mesh!` call per surface would give a detector a few hundred plot objects; concatenating shifts
that to one draw call per translucency class, which is what keeps CairoMakie's per-primitive sorting
tractable (plan §6).
"""
struct MeshBuffer
    points::Vector{Makie.Point3f}
    faces::Vector{Makie.GeometryBasics.GLTriangleFace}
end

MeshBuffer() = MeshBuffer(Makie.Point3f[], Makie.GeometryBasics.GLTriangleFace[])

function Base.append!(buf::MeshBuffer, m::Mesh)
    offset = length(buf.points)
    append!(buf.points, mesh_points(m))
    for f in triangle_faces(m.connections)
        push!(buf.faces, Makie.GeometryBasics.GLTriangleFace(f[1] + offset, f[2] + offset, f[3] + offset))
    end
    buf
end

Base.isempty(buf::MeshBuffer) = isempty(buf.faces)

geometry_mesh(buf::MeshBuffer) = Makie.GeometryBasics.Mesh(buf.points, buf.faces)

"""
    line_segments(m::CSG.Mesh) -> Vector{Point3f}

Flatten a mesh whose `connections` are index *pairs* into the point vector `linesegments!` expects.

Use it on `mesh(s, n_arc, n_vert_lines)`; on a face mesh use [`face_outlines`](@ref) instead.
"""
function line_segments(m::Mesh)
    pts = mesh_points(m)
    segs = Makie.Point3f[]
    for c in m.connections, i in 1:(length(c) - 1)
        push!(segs, pts[c[i]], pts[c[i + 1]])
    end
    segs
end

"""
    face_outlines(m::CSG.Mesh) -> Vector{Point3f}

Flatten the closed outline of every face of a mesh into `linesegments!` point pairs.

This is the grid `seriestype = :mesh3d` draws over its faces. Drawing it explicitly rather than with
`wireframe!` is what keeps the fan triangulation invisible: `wireframe!` would also draw the
diagonals [`triangle_faces`](@ref) introduces.
"""
function face_outlines(m::Mesh)
    pts = mesh_points(m)
    segs = Makie.Point3f[]
    for f in m.connections
        length(f) < 2 && continue
        for i in eachindex(f)
            push!(segs, pts[f[i]], pts[f[i % length(f) + 1]])
        end
    end
    segs
end

"""
    point_coordinates(pts) -> Vector{Point3f}

Convert CSG points of either coordinate system into display-unit Makie points.
"""
point_coordinates(pts::AbstractVector{<:CylindricalPoint}) =
    point_coordinates([CartesianPoint(p) for p in pts])

function point_coordinates(pts::AbstractVector{<:CartesianPoint})
    f = display_factor(internal_length_unit)
    Makie.Point3f[Makie.Point3f(f * p.x, f * p.y, f * p.z) for p in pts]
end

point_coordinates(p::CartesianPoint) = point_coordinates([p])
point_coordinates(p::CylindricalPoint) = point_coordinates([CartesianPoint(p)])
