# Backend-agnostic plotting helpers of the CSG submodule.
#
# The counterpart of `src/PlotRecipes/plot_utils.jl` one module up: data preparation and label text
# that both the Plots.jl recipes under `plotting/` and the Makie extension need. They live *outside*
# `plotting/` because that whole directory is deleted when the Plots recipes go (migration plan
# step 11), while these have to survive it.

## Legend labels
#
# Plots builds its legend from every series' `label`; a Makie recipe cannot produce legend entries at
# all, so its figure-level wrapper builds them from explicit elements. Both must print the same text,
# which is why the strings live here and not in either recipe. Surface primitives carry their own
# `get_label_name` methods next to their definitions (`"Cone Mantle"`, `"Polygon"`, ...).

"""
    get_label_name(geometry) -> String

Give the legend label of a geometry.

A volume primitive is labelled by its type name (`"Box"`, `"Cone"`), a surface primitive by its
descriptive name (`"Cone Mantle"`), and any union, intersection or difference by `"CSG"`.
"""
get_label_name(p::AbstractVolumePrimitive) = string(nameof(typeof(p)))
get_label_name(::AbstractConstructiveGeometry) = "CSG"
get_label_name(::Ellipse) = "Ellipse"
get_label_name(::Edge) = "Edge"


## Cross sections

"""
    get_crosssection_and_value(x, y, z, φ) -> (axis, value)

Resolve the cross-section plotting keywords into the cut axis and its position in internal units.

Exactly one of the four is expected to be given; if none is, the cut defaults to `y = 0`.
"""
function get_crosssection_and_value(x::Union{Real, LengthQuantity, Missing}, y::Union{Real, LengthQuantity, Missing}, z::Union{Real, LengthQuantity, Missing}, φ::Union{Real, AngleQuantity, Missing})::Tuple{Symbol, Real}
    crosssections = Dict(:x => x, :y => y, :z => z, :φ => φ)
    idx = findfirst(cs -> !ismissing(cs), crosssections)
    if isnothing(idx)
        :y, 0
    else
        return idx, to_internal_units(crosssections[idx])
    end
end

"""
    is_coplanar(s::AbstractSurfacePrimitive, plane::Tuple{Symbol, T}) -> Bool

Test whether the surface lies in the cutting plane rather than crossing it.

A cross section through such a surface is sampled twice as finely, because the samples landing on
the plane are then the whole picture.
"""
function is_coplanar(s::AbstractSurfacePrimitive{T}, plane::Tuple{Symbol,T}) where {T}
    axis, slice = plane
    crosssections = Dict(   :x => CartesianVector{T}(1,0,0),
                            :y => CartesianVector{T}(0,1,0),
                            :z => CartesianVector{T}(0,0,1),
                            :φ => CartesianVector{T}(cos(π/2+slice),sin(π/2+slice),0)
                   )
    pointtype = axis == :φ ? CylindricalPoint : CartesianPoint
    point = pointtype(vertices(s,1)[1])
    abs(abs(dot(normal(s), crosssections[axis])) - 1) < csg_default_tol(T) && abs(getproperty(point, axis) - slice) < csg_default_tol(T)
end


## Sampling

"""
    sample(p::AbstractVolumePrimitive, spacing) -> Vector{CartesianPoint}

Sample the surface of a volume primitive at points spaced roughly `spacing` apart.

This is the marker-based representation of a geometry (`:samplesurface`, `:slice`), and the only one
that stays exact for differences and intersections, where the mesh-based ones draw whole primitives.
"""
function sample(p::AbstractVolumePrimitive{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    fs = surfaces(p)
    vs = Vector{CartesianPoint{T}}()
    for s in fs
        append!(vs, sample(s, spacing))
    end
    vs
end
