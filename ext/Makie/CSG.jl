# Makie recipe and figure-level wrapper for the constructive solid geometry layer: volume primitives,
# surface primitives, their unions, intersections and differences, plus the points, vectors, edges
# and ellipses that go with them. 31 of the ~142 parity cases, and the mirror image of the recipes
# under `src/ConstructiveSolidGeometry/plotting/`.
#
# Two layers as everywhere else (plan §1.1 / D1): `CSGPlot` draws into whatever axis it is handed,
# `csgplot(geom; …)` builds the `Figure`, the `Axis3`/`Axis` and the `Legend`.
#
# The `seriestype` keyword is kept, spelling and all. It is Plots vocabulary and Makie has no such
# attribute, but it is what selects the representation at all 21 doc call sites and in the case
# matrix, and there is nothing to gain from renaming it mid-migration.
#
# What is *not* a literal translation:
#
#   * **`:mesh3d` + `connections` becomes two idiomatic calls.** A face list goes to `mesh!` through
#     `GeometryBasics.Mesh`, a line list to `linesegments!` (plan §1.8). GR's painter-ordered
#     `mesh3d` is replaced by a real, depth-tested surface, which is one of the places Makie output
#     should be visibly better.
#   * **One draw call per translucency class**, not one per surface. The Plots recipes emit two
#     `@series` per face; here every closed face of a geometry lands in one mesh and every subtracted
#     one in another, because a detector otherwise reaches a few hundred plot objects.
#   * **Backend sniffing is gone.** The six `occursin("GRBackend", …)` guards existed to work around
#     GR-vs-PyPlot differences; the equal aspect they gated is applied unconditionally here.

"""
    csgplot(geometry; kwargs...)
    csgplot!(ax, geometry; kwargs...)

Plot a CSG geometry: a volume or surface primitive, a union, intersection or difference of them, an
`Edge`, an `Ellipse`, a `CartesianPoint`, a vector of points, or a `CartesianVector`.

The non-mutating form returns a `FigureAxisPlot` carrying the axis labels and the legend; the
mutating form draws into an existing axis and sets neither. `plot(geometry)` routes here.

## Keywords
* `seriestype`: the representation. `:csg` (default) draws the faces with a wireframe on top,
  `:mesh3d` the faces with their cell grid, `:wireframe` the wireframe alone, `:samplesurface`
  markers sampled over the surface, and `:slice` a two-dimensional cut.
* `x`, `y`, `z`, `φ`: with `seriestype = :slice`, the position of the cut; the default is `y = 0`.
  A `φ` cut shows only `r >= 0` unless `full_det = true`.
* `n_arc`: number of points per full arc of a curved surface, i.e. the mesh resolution.
* `n_samples`: number of samples per `CSG_scale` for `:samplesurface` and `:slice`. Raising it is
  what turns a sparse marker cloud into a solid-looking surface, at proportional cost.
* `color`, `alpha`, `linecolor`, `linewidth`, `markersize`: drawing attributes, each defaulting to
  the value for the chosen representation.

Differences and intersections are exact only under `:samplesurface` and `:slice`. The mesh-based
representations draw whole primitives, and show a subtracted one as a thin white shell instead.
"""
@recipe CSGPlot (geometry,) begin
    "Representation: `:csg`, `:mesh3d`, `:wireframe`, `:samplesurface` or `:slice`."
    seriestype = :csg
    "Colour of the faces and markers; `Makie.automatic` takes the first palette entry."
    color = Makie.automatic
    "Opacity; `Makie.automatic` takes the value for the chosen representation."
    alpha = Makie.automatic
    "Colour of the wireframe and of the `:mesh3d` cell grid; `Makie.automatic` is black."
    linecolor = Makie.automatic
    "Width of the wireframe lines."
    linewidth = Makie.automatic
    "Size of the markers of `:samplesurface` and `:slice`."
    markersize = Makie.automatic
    "Points per full arc of a curved surface."
    n_arc = 40
    "Number of lines drawn along a curved surface by the wireframe."
    n_vert_lines = 2
    "Samples per `CSG_scale`; `Makie.automatic` is 40 in 3D and 100 or 200 for a slice."
    n_samples = Makie.automatic
    "Length scale the sample spacing is derived from; `missing` takes the geometry's own extent."
    CSG_scale = missing
    "Position of a `:slice` cut."
    x = missing
    "Position of a `:slice` cut."
    y = missing
    "Position of a `:slice` cut."
    z = missing
    "Angle of a `:slice` cut; the half plane `r >= 0` unless `full_det = true`."
    φ = missing
    "Draw both halves of a `φ` slice."
    full_det = false
    "`:polar` draws a `z` slice in (θ, r), to overlay a `PolarAxis`; `:none` draws it in (x, y)."
    projection = :none
    "Widen a dimension the geometry has no extent in, so equal data scaling cannot collapse it."
    pad_flat = true
end

# Everything drawn as a surface: the two primitive families and any combination of them.
const CSGGeometry = Union{AbstractVolumePrimitive, AbstractSurfacePrimitive,
                          AbstractConstructiveGeometry}
const CSGPointLike = Union{AbstractCoordinatePoint, AbstractVector{<:AbstractCoordinatePoint}}

# `plot(obj)` / `plot!(ax, obj)` route to the recipe, which is what keeps the doc call sites working
# (D4). `convert_arguments` is the identity: the recipe consumes the geometry itself.
for GT in (:CSGGeometry, :CSGPointLike, :AbstractLinePrimitive, :CartesianVector,
           :(AbstractVector{<:AbstractSurfacePrimitive}), :(AbstractVector{<:AbstractLinePrimitive}))
    @eval Makie.plottype(::$GT) = CSGPlot
    @eval Makie.convert_arguments(::Type{<:CSGPlot}, g::$GT) = (g,)
end

# A slice is a plane cut and lives on a plain `Axis`; every other representation is geometry in
# space. `preferred_axis_type` runs before the recipe's defaults are merged, so the raw kwargs are
# the only place `seriestype` exists yet.
_representation(p) = Makie.to_value(get(p.kw, :seriestype, :csg))
Makie.preferred_axis_type(p::CSGPlot{<:Tuple{<:CSGGeometry}}) =
    _representation(p) === :slice ? Makie.Axis : Makie.Axis3
Makie.preferred_axis_type(::CSGPlot) = Makie.Axis3


## ----------------------------------------------------------------------------------------------
## Decomposition
## ----------------------------------------------------------------------------------------------

# The precision type of a geometry. Not `eltype`: a primitive is not an array, so `eltype` gives
# `Any` and every `T(…)` conversion below would fail.
_csg_eltype(::AbstractGeometry{T}) where {T} = T
_csg_eltype(v::AbstractVector{<:AbstractGeometry}) = _csg_eltype(first(v))

# The primitives a geometry is drawn as, each with the flag that says whether it is *added* to the
# result or *removed* from it. A removed primitive is drawn as a thin white shell, which is the whole
# reason the mesh-based representations can show a difference at all.
_csg_parts(p::AbstractPrimitive) = [(p, true)]
_csg_parts(p::AbstractVolumePrimitive) = [(p, isClosedPrimitive(p))]
_csg_parts(csg::AbstractConstructiveGeometry) = [(q, isClosedPrimitive(q)) for q in primitives(csg)]

# The surfaces of a primitive: its own faces, or itself when it already is a face.
_csg_faces(p::AbstractVolumePrimitive) = collect(surfaces(p))
_csg_faces(s::AbstractSurfacePrimitive) = [s]

# The scale the sample spacing is measured in. `get_scale` is `extremum` for a primitive and the
# extent of the whole tree for a combination, so a subtracted cone is sampled as finely as the box
# it is cut out of.
_csg_scale(p::AbstractPrimitive) = extremum(p)
_csg_scale(csg::AbstractConstructiveGeometry) = get_scale(csg)

# Samples on the surface of one primitive that survive the geometry: a difference removes the part of
# the box's surface that lies inside the cone, and only `in` knows that.
_csg_filter(::AbstractPrimitive) = _ -> true
_csg_filter(csg::AbstractConstructiveGeometry) = pt -> pt in csg


## ----------------------------------------------------------------------------------------------
## The recipe
## ----------------------------------------------------------------------------------------------

function Makie.plot!(p::CSGPlot{<:Tuple{<:CSGGeometry}})
    st = p.seriestype[]
    if st === :slice
        _plot_csg_slice!(p)
    elseif st === :samplesurface
        _plot_csg_samples!(p)
    elseif st in (:csg, :mesh3d, :wireframe)
        _plot_csg_mesh!(p, st)
    else
        throw(ArgumentError("unknown seriestype $(repr(st)); use :csg, :mesh3d, :wireframe, " *
                            ":samplesurface or :slice"))
    end
    st !== :slice && p.pad_flat[] && pad_flat_extent!(p)
    return p
end

_face_color(p) = p.color[] === Makie.automatic ? contact_color(1) : p.color[]

function _plot_csg_mesh!(p, st::Symbol)
    geom = p[1][]
    n_arc, n_vert = p.n_arc[], p.n_vert_lines[]
    solid, shell = MeshBuffer(), MeshBuffer()
    wires, thin_wires = Makie.Point3f[], Makie.Point3f[]
    for (prim, closed) in _csg_parts(geom), s in _csg_faces(prim)
        # `:wireframe` never needs the face mesh, and the two line representations come from
        # different meshes, so exactly one `mesh(…)` call happens per surface per representation.
        faces = st === :wireframe ? nothing : mesh(s, n_arc)
        faces === nothing || append!(closed ? solid : shell, faces)
        segments = st === :mesh3d ? face_outlines(faces) : line_segments(mesh(s, n_arc, n_vert))
        append!(closed ? wires : thin_wires, segments)
    end

    # A subtracted primitive is a white shell so the solid it was cut from stays readable through it;
    # `fillcolor := :white` in the Plots recipe, and the same idea.
    facealpha = p.alpha[] === Makie.automatic ? SSD_ALPHAS.mesh : p.alpha[]
    isempty(solid) || Makie.mesh!(p, geometry_mesh(solid);
                                  color = _face_color(p), alpha = facealpha, shading = Makie.NoShading)
    isempty(shell) || Makie.mesh!(p, geometry_mesh(shell);
                                  color = :white, alpha = SSD_ALPHAS.subtracted,
                                  shading = Makie.NoShading)

    linecolor = p.linecolor[] === Makie.automatic ?
        :black : p.linecolor[]
    lw = p.linewidth[] === Makie.automatic ?
        (st === :mesh3d ? SSD_LINEWIDTHS.meshgrid : SSD_LINEWIDTHS.wireframe) : p.linewidth[]
    isempty(wires) || Makie.linesegments!(p, wires; color = linecolor, linewidth = lw)
    isempty(thin_wires) ||
        Makie.linesegments!(p, thin_wires; color = linecolor, linewidth = SSD_LINEWIDTHS.subtracted)
    p
end

function _plot_csg_samples!(p)
    geom = p[1][]
    n = p.n_samples[] === Makie.automatic ? 40 : p.n_samples[]
    spacing = _sample_spacing(p, geom, n)
    keep = _csg_filter(geom)
    pts = CartesianPoint{_csg_eltype(geom)}[]
    for (prim, _) in _csg_parts(geom)
        append!(pts, filter(keep, sample(prim, spacing)))
    end
    Makie.scatter!(p, point_coordinates(pts);
                   color = _face_color(p),
                   alpha = p.alpha[] === Makie.automatic ? SSD_ALPHAS.samplesurface : p.alpha[],
                   markersize = _markersize(p, SSD_MARKERSIZES.samplesurface),
                   strokewidth = 0)
    p
end

_markersize(p, default) = p.markersize[] === Makie.automatic ? default : p.markersize[]

function _sample_spacing(p, geom, n_samples)
    T = _csg_eltype(geom)
    scale = ismissing(p.CSG_scale[]) ? _csg_scale(geom) : p.CSG_scale[]
    T(scale / n_samples)
end

# The number of samples a slice defaults to. A surface lying *in* the cutting plane contributes its
# whole area to the picture rather than a single line of points, so half the density suffices there.
function _slice_n_samples(geom, axis::Symbol, slice)
    T = _csg_eltype(geom)
    planar = filter(s -> s isa AbstractPlanarSurfacePrimitive{T}, collect(surfaces(geom)))
    any(s -> is_coplanar(s, (axis, T(slice))), planar) ? 100 : 200
end

# The two Cartesian coordinates a cut is drawn in, and the names they carry. A `φ` cut is taken after
# rotating the geometry onto the x-z half plane, so its coordinates are `x` and `z` while its axes
# are labelled `r` and `z`.
function _slice_projection(axis::Symbol)
    pointtype = axis === :φ ? CylindricalPoint : CartesianPoint
    keep = findall(!=(axis), collect(fieldnames(pointtype)))
    coords = collect(fieldnames(CartesianPoint))[keep]
    names = collect(fieldnames(pointtype))[keep]
    (coords[1], coords[2]), (String(names[1]), String(names[2]))
end

function _plot_csg_slice!(p)
    geom = p[1][]
    T = _csg_eltype(geom)
    axis, slice = get_crosssection_and_value(p.x[], p.y[], p.z[], p.φ[])
    n = p.n_samples[] === Makie.automatic ? _slice_n_samples(geom, axis, slice) : p.n_samples[]

    # A `φ` cut is a rotation away from a `y = 0` cut, which is how both backends do it: rotate the
    # geometry by -φ and take the plane `y = 0` through it. The sample spacing is halved because the
    # cut then has to catch points on a plane rather than in a slab.
    axisrot, position = axis, T(slice)
    if axis === :φ
        geom = _rotate_to_φ(geom, T(slice))
        axisrot, position = :y, zero(T)
    end
    spacing = _sample_spacing(p, geom, n)
    axis === :φ && (spacing *= T(0.5))

    keep = _csg_filter(geom)
    on_plane(pt) = abs(getproperty(pt, axisrot) - position) < spacing / 2
    pts = CartesianPoint{T}[]
    for (prim, _) in _csg_parts(geom)
        append!(pts, filter(pt -> keep(pt) && on_plane(pt), sample(prim, spacing)))
    end

    f = display_factor(internal_length_unit)
    if p.projection[] === :polar
        # A `PolarAxis` takes (θ, r), with θ in radians (D7). Only a `z` cut has a polar form: the
        # other cuts are planes through the axis, not discs around it.
        axis === :z || throw(ArgumentError(
            "projection = :polar is only defined for a `z` slice, not a `$(axis)` one"))
        cyl = CylindricalPoint.(pts)
        us = [pt.φ for pt in cyl]
        vs = [f * pt.r for pt in cyl]
        return _draw_slice_markers!(p, us, vs)
    end
    (cu, cv), _ = _slice_projection(axis)
    us = [f * getproperty(pt, cu) for pt in pts]
    vs = [f * getproperty(pt, cv) for pt in pts]
    # Only the half plane the cut was asked for, unless both are wanted: the rotated geometry extends
    # to negative `x`, which is the φ + 180° cross section.
    if axis === :φ && !p.full_det[]
        keep_idx = findall(>=(0), us)
        us, vs = us[keep_idx], vs[keep_idx]
    end

    _draw_slice_markers!(p, us, vs)
end

# The markers of a slice, in whichever pair of coordinates the projection produced.
function _draw_slice_markers!(p, us, vs)
    Makie.scatter!(p, us, vs;
                   color = _face_color(p),
                   alpha = p.alpha[] === Makie.automatic ? 1.0 : p.alpha[],
                   markersize = _markersize(p, SSD_MARKERSIZES.slice),
                   strokewidth = 0)
    p
end

_rotate_to_φ(p::AbstractPrimitive, φ) = rotate(p, typeof(φ).(RotZ(-φ)))
_rotate_to_φ(csg::AbstractConstructiveGeometry{T}, φ) where {T} =
    transform(csg, (rotation = T.(RotZ(-φ)), translation = CartesianVector{T}()))


## ----------------------------------------------------------------------------------------------
## Lines, points and vectors
## ----------------------------------------------------------------------------------------------

# An `Ellipse` is drawn as the `n_arc` edges it decomposes into, an `Edge` as itself. Plots emitted
# one series per edge, cycling the palette while forcing every line black; one `linesegments!` call
# draws the same picture.
function Makie.plot!(p::CSGPlot{<:Tuple{<:AbstractLinePrimitive}})
    _plot_line_primitives!(p, _line_edges(p[1][], p.n_arc[]))
end

Makie.plot!(p::CSGPlot{<:Tuple{<:AbstractVector{<:AbstractLinePrimitive}}}) =
    _plot_line_primitives!(p, reduce(vcat, (_line_edges(l, p.n_arc[]) for l in p[1][]); init = Edge[]))

_line_edges(e::Edge, ::Int) = [e]
_line_edges(e::Ellipse, n_arc::Int) = collect(edges(e, n = n_arc))

function _plot_line_primitives!(p, es)
    pts = Makie.Point3f[]
    for e in es
        append!(pts, point_coordinates([e.a, e.b]))
    end
    Makie.linesegments!(p, pts;
                        color = p.linecolor[] === Makie.automatic ? :black : p.linecolor[],
                        linewidth = p.linewidth[] === Makie.automatic ?
                            SSD_LINEWIDTHS.wireframe : p.linewidth[])
    p.pad_flat[] && pad_flat_extent!(p)
    p
end

# A vector of surfaces is what `surfaces(primitive)` returns; drawing it goes through the same mesh
# path as the primitive it came from.
function Makie.plot!(p::CSGPlot{<:Tuple{<:AbstractVector{<:AbstractSurfacePrimitive}}})
    st = p.seriestype[]
    st === :slice && throw(ArgumentError(":slice is not defined for a vector of surfaces"))
    st === :samplesurface ? _plot_csg_samples!(p) : _plot_csg_mesh!(p, st)
    p.pad_flat[] && pad_flat_extent!(p)
    p
end

_csg_parts(v::AbstractVector{<:AbstractSurfacePrimitive}) = [(s, true) for s in v]
_csg_scale(v::AbstractVector{<:AbstractSurfacePrimitive}) = maximum(extremum, v)

# The adaptive marker size of the Plots recipe, kept because it is backend-agnostic: a cloud of a few
# thousand points needs smaller markers than a handful, or it is a solid blob.
function _auto_markersize(n::Integer)
    default, n_thresh = 4.0, 10
    invscale = sqrt(max(1, 2 + (n - 1 - n_thresh) / n_thresh))
    SSD_MARKERSIZES.point_scale * max(default / invscale, 0.1)
end

function Makie.plot!(p::CSGPlot{<:Tuple{<:CSGPointLike}})
    pts = point_coordinates(_as_point_vector(p[1][]))
    Makie.scatter!(p, pts;
                   color = _face_color(p),
                   markersize = _markersize(p, _auto_markersize(length(pts))),
                   alpha = p.alpha[] === Makie.automatic ? 1.0 : p.alpha[])
    p.pad_flat[] && pad_flat_extent!(p)
    p
end

_as_point_vector(pt::AbstractCoordinatePoint) = [pt]
_as_point_vector(v::AbstractVector{<:AbstractCoordinatePoint}) = v

# `plot(::CartesianVector)` never worked under Plots: its `:vector` series type runs after the
# positional arguments have been converted to series data, and `(CartesianPoint, CartesianVector)`
# has no such conversion, so the hand-built 13-point arrowhead was unreachable (plan F2). Makie
# draws an arrow properly, which also retires the arrowhead construction.
function Makie.plot!(p::CSGPlot{<:Tuple{<:CartesianVector}})
    v = p[1][]
    f = display_factor(internal_length_unit)
    Makie.arrows3d!(p, [Makie.Point3f(0, 0, 0)], [Makie.Vec3f(f * v.x, f * v.y, f * v.z)];
                    color = p.linecolor[] === Makie.automatic ? :black : p.linecolor[],
                    shaftradius = 0.012, tipradius = 0.04, tiplength = 0.12,
                    shading = Makie.NoShading)
    p
end


## ----------------------------------------------------------------------------------------------
## Figure-level wrapper
## ----------------------------------------------------------------------------------------------

"""
    geometry_axis_setup(; seriestype, x, y, z, φ, kwargs...) -> NamedTuple

Give the `ssd_axis` keywords for a plot of geometry in space.

Depends on nothing but the representation and the cut, so the detector wrapper shares it with the
CSG one: what a geometry plot needs from its axis is the same either way.
"""
function geometry_axis_setup(; seriestype = :csg, x = missing, y = missing, z = missing,
                             φ = missing, kwargs...)
    lu = internal_length_unit
    if seriestype === :slice
        axis, _ = get_crosssection_and_value(x, y, z, φ)
        _, (xname, yname) = _slice_projection(axis)
        (; axis = (; axistype = Makie.Axis, xname, yname, xunit = lu, yunit = lu,
                   aspect_ratio_one = true))
    else
        # `aspect = :data` on every three-dimensional plot: a torus of radius 2 with a 0.5 tube is
        # flat, and GR drew it in a cube whatever its proportions. Accepted as an improvement under
        # plan A1, and it is what the Plots point and line recipes asked for anyway.
        (; axis = (; axistype = Makie.Axis3, xname = "x", yname = "y", zname = "z",
                   xunit = lu, yunit = lu, zunit = lu, aspect_ratio_one = true))
    end
end

# The legend entry Plots produced from the label of the first series, as an explicit element (§1.3).
# Points carry no entry: their Plots label is the auto-generated `"y1"`, which is not a legend entry
# on either side.
function _csg_legend_entries(geom; seriestype = :csg, color = Makie.automatic, kwargs...)
    facecolor = color === Makie.automatic ? contact_color(1) : color
    # A bare geometry's wireframe is black, unlike a detector part's, which takes the part's colour.
    swatch = seriestype === :wireframe ? :black : facecolor
    (Makie.LegendElement[representation_element(seriestype, swatch)],
     String[get_label_name(geom)])
end

_csg_legend_entries(::CSGPointLike; kwargs...) = (Makie.LegendElement[], String[])
_csg_legend_entries(l::AbstractLinePrimitive; linecolor = Makie.automatic, kwargs...) =
    (Makie.LegendElement[representation_element(:wireframe, color_or_default(linecolor, :black))],
     String[get_label_name(l)])
_csg_legend_entries(::CartesianVector; kwargs...) = (Makie.LegendElement[], String[])
_csg_legend_entries(v::AbstractVector{<:AbstractSurfacePrimitive}; kwargs...) =
    (Makie.LegendElement[geometry_legend_element(contact_color(1))], String["Faces"])

const CSGPlottable = Union{CSGGeometry, CSGPointLike, AbstractLinePrimitive, CartesianVector,
                           AbstractVector{<:AbstractSurfacePrimitive},
                           AbstractVector{<:AbstractLinePrimitive}}

# A method of the recipe's own function, more specific than the `(args...)` method `@recipe`
# generates, so this is what a caller of `csgplot(geom)` gets.
function csgplot(gp::Makie.GridPosition, geom::CSGPlottable; axis = (;), legend = (;), kwargs...)
    setup = geometry_axis_setup(; kwargs...)
    gl = panel_layout(gp)
    ax = ssd_axis(gl[1, 1]; setup.axis..., axis...)
    plt = csgplot!(ax, geom; kwargs...)
    finish_geometry_axis!(ax, plt, axis)
    legend === false || add_legend!(gl, _csg_legend_entries(geom; kwargs...)...;
                                    corner = legend_corner(ax, plt), legend...)
    Makie.AxisPlot(ax, plt)
end

function csgplot(geom::CSGPlottable; figure = (;), kwargs...)
    fig = Makie.Figure(; figure...)
    ap = csgplot(fig[1, 1], geom; kwargs...)
    Makie.FigureAxisPlot(fig, ap.axis, ap.plot)
end

# `plot(geom)` keeps producing the labelled figure with a legend that Plots produced, rather than a
# bare mesh; `plot!(ax, geom)` stays on the recipe, since an existing axis belongs to its caller.
Makie.plot(geom::CSGPlottable; kwargs...) = csgplot(geom; kwargs...)
