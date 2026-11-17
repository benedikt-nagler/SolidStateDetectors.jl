# Makie recipe and figure-level wrapper for a detector and its parts: `SolidStateDetector`,
# `Contact`, `Semiconductor` and `Passive`. 18 of the ~142 parity cases.
#
# This layer owns almost no geometry code. A detector *is* a list of parts, each of which is a CSG
# geometry with a colour and an opacity, so the recipe delegates every part to `csgplot!` (step 5)
# and contributes exactly three things the CSG layer cannot know about:
#
#   * **which parts** are drawn and in what order (`detector_parts`), shared with the legend so the
#     two can never disagree;
#   * **the colours and opacities** of the parts, which is what makes a detector readable: contacts
#     and passives at `alpha = 0.2` over a semiconductor at `0.5`, each in its own colour, with the
#     wireframe taking the part's colour rather than the black a bare geometry gets;
#   * **one shared sampling scale.** The marker-based representations divide `CSG_scale` by
#     `n_samples`, so sampling each part against its own extent would draw a 3 mm point contact at a
#     tenth of the spacing of the 40 mm crystal, and the contact would come out as a solid blob. The
#     Plots recipe takes `maximum(get_scale(...))` over all parts; so does this one.
#
# The legend is where the two backends differ most (plan §1.3). Every Plots recipe here emits an
# **empty `:path` series** per part, purely so the legend gets an entry drawn at `alpha = 1` while
# the geometry itself is drawn at `0.2`. Makie collects entries only from direct children of the
# axis, so a recipe that draws all parts inside itself produces none; the wrapper builds the `Legend`
# from explicit elements instead, and the empty-series hack disappears rather than being ported.

"""
    detectorplot(detector; kwargs...)
    detectorplot!(ax, detector; kwargs...)

Plot a [`SolidStateDetector`](@ref) or one of its parts: a `Contact`, the `Semiconductor`, or a
`Passive`.

The non-mutating form returns a `FigureAxisPlot` carrying the axis labels and the legend; the
mutating form draws into an existing axis and sets neither. `plot(detector)` routes here.

## Keywords
* `seriestype`: the representation, as for [`csgplot`](@ref). `:csg` (default), `:mesh3d`,
  `:wireframe`, `:samplesurface` or `:slice`.
* `show_semiconductor`: draw the semiconductor as well as the contacts. Off by default, because it
  encloses everything else.
* `show_passives`: draw the passives. On by default.
* `x`, `y`, `z`, `φ`: with `seriestype = :slice`, the position of the cut.
* `n_arc`, `n_samples`, `linewidth`, `markersize`: as for [`csgplot`](@ref).

Contacts sit inside the semiconductor, so `show_semiconductor = true` is only readable under
`:wireframe` or a `:slice`; the mesh representations bury the contacts behind it.
"""
@recipe DetectorPlot (detector,) begin
    "Representation: `:csg`, `:mesh3d`, `:wireframe`, `:samplesurface` or `:slice`."
    seriestype = :csg
    "Draw the semiconductor, which encloses the contacts."
    show_semiconductor = false
    "Draw the passives."
    show_passives = true
    "Points per full arc of a curved surface."
    n_arc = 40
    "Number of lines drawn along a curved surface by the wireframe."
    n_vert_lines = 2
    "Samples per `CSG_scale`; `Makie.automatic` is 40 in 3D and 100 or 200 for a slice."
    n_samples = Makie.automatic
    "Width of the wireframe lines."
    linewidth = Makie.automatic
    "Size of the markers of `:samplesurface` and `:slice`."
    markersize = Makie.automatic
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
    "Colour of every part; `automatic` gives each its own (contacts by `id`)."
    color = Makie.automatic
    "Opacity of every part; `automatic` keeps the per-part default."
    alpha = Makie.automatic
    "`:polar` draws a `z` slice in (θ, r), to overlay a potential drawn on a `PolarAxis`."
    projection = :none
end

const DetectorPart = Union{Contact, Semiconductor, Passive}

for GT in (:SolidStateDetector, :DetectorPart)
    @eval Makie.plottype(::$GT) = DetectorPlot
    @eval Makie.convert_arguments(::Type{<:DetectorPlot}, d::$GT) = (d,)
end

Makie.preferred_axis_type(p::DetectorPlot) =
    Makie.to_value(get(p.kw, :seriestype, :csg)) === :slice ? Makie.Axis : Makie.Axis3


## ----------------------------------------------------------------------------------------------
## The parts of a detector
## ----------------------------------------------------------------------------------------------

"""
    detector_parts(det::SolidStateDetector; show_semiconductor = false, show_passives = true)

List the parts of a detector in the order they are drawn and listed in the legend: contacts first,
then the semiconductor if shown, then the passives.

One definition of the order, consumed by both the recipe and `legend_entries`.
"""
function detector_parts(det::SolidStateDetector;
                        show_semiconductor::Bool = false, show_passives::Bool = true)
    parts = DetectorPart[]
    append!(parts, det.contacts)
    show_semiconductor && push!(parts, det.semiconductor)
    show_passives && !ismissing(det.passives) && append!(parts, det.passives)
    parts
end

detector_parts(part::DetectorPart; kwargs...) = DetectorPart[part]

# Colour and opacity per part. Both are pinned rather than cycled (D3): `Cycled(i)` would index
# Makie's Wong palette, where the Plots recipes index the palette `SSD_PALETTE` reproduces.
part_color(c::Contact) = contact_color(c.id)
part_color(::Semiconductor) = SSD_COLORS.semiconductor
part_color(::Passive) = SSD_COLORS.passive

part_alpha(::Contact) = SSD_ALPHAS.contact
part_alpha(::Semiconductor) = SSD_ALPHAS.semiconductor
part_alpha(::Passive) = SSD_ALPHAS.passive


## ----------------------------------------------------------------------------------------------
## The recipe
## ----------------------------------------------------------------------------------------------

function Makie.plot!(p::DetectorPlot)
    parts = detector_parts(p[1][];
                           show_semiconductor = p.show_semiconductor[],
                           show_passives = p.show_passives[])
    st = p.seriestype[]
    # One spacing for the whole detector, or a 3 mm contact is sampled ten times finer than the
    # crystal around it and reads as a solid blob next to a sparse cloud.
    scale = st in (:samplesurface, :slice) ?
        maximum(get_scale(part.geometry) for part in parts) : missing
    for part in parts
        # An explicit colour or opacity overrides the per-part scheme for the whole detector, which
        # is what makes a single-colour slice or a faded semiconductor possible.
        color = color_or_default(p.color[], part_color(part))
        # A slice is an outline, not a body: drawn at the part's 0.2 it all but disappears.
        alpha = st === :slice ? p.alpha[] : color_or_default(p.alpha[], part_alpha(part))
        csgplot!(p, part.geometry;
                 seriestype = st, color, alpha,
                 # A detector part draws its wireframe in its own colour, where a bare geometry
                 # draws it black: that is what tells two contacts apart in a `:wireframe` plot.
                 linecolor = color, CSG_scale = scale,
                 n_arc = p.n_arc[], n_vert_lines = p.n_vert_lines[], n_samples = p.n_samples[],
                 linewidth = p.linewidth[], markersize = p.markersize[],
                 x = p.x[], y = p.y[], z = p.z[], φ = p.φ[], full_det = p.full_det[],
                 projection = p.projection[], pad_flat = false)
    end
    st === :slice || pad_flat_extent!(p)
    return p
end


## ----------------------------------------------------------------------------------------------
## Figure-level wrapper
## ----------------------------------------------------------------------------------------------

const DetectorPlottable = Union{SolidStateDetector, DetectorPart}

function detectorplot(gp::Makie.GridPosition, det::DetectorPlottable;
                      axis = (;), legend = (;), kwargs...)
    # The axis of a detector is the axis of the geometry it is made of, from the same helper the CSG
    # wrapper uses.
    setup = geometry_axis_setup(; kwargs...)
    gl = panel_layout(gp)
    ax = ssd_axis(gl[1, 1]; setup.axis..., axis...)
    plt = detectorplot!(ax, det; kwargs...)
    finish_geometry_axis!(ax, plt, axis)
    st = get(kwargs, :seriestype, :csg)
    parts = detector_parts(det;
                           show_semiconductor = get(kwargs, :show_semiconductor, false),
                           show_passives = get(kwargs, :show_passives, true))
    given = get(kwargs, :color, Makie.automatic)
    legend === false ||
        add_legend!(gl,
                    Makie.LegendElement[representation_element(st, color_or_default(given, part_color(q)))
                                        for q in parts],
                    String[plot_label(q) for q in parts]; corner = legend_corner(ax, plt), legend...)
    Makie.AxisPlot(ax, plt)
end

function detectorplot(det::DetectorPlottable; figure = (;), kwargs...)
    fig = Makie.Figure(; figure...)
    ap = detectorplot(fig[1, 1], det; kwargs...)
    Makie.FigureAxisPlot(fig, ap.axis, ap.plot)
end

Makie.plot(det::DetectorPlottable; kwargs...) = detectorplot(det; kwargs...)
