# Makie recipe and figure-level wrapper for charge drift paths: an `EHDriftPath` or a vector of them.
# 3 of the ~142 parity cases.
#
# The picture is four strokes: the electron trajectories in red, the hole trajectories in green, and
# a marker at the end of each, where the carrier was collected. The Plots recipe emits those four
# series *per path* and a vector of paths therefore produces 4n series; here each species is drawn in
# one `lines!` with a NaN vertex between paths, so an event with a thousand carriers stays four plot
# objects. Makie breaks a line at a non-finite point, which is what makes that legitimate rather than
# a trick, and the same "one draw call per class" choice the CSG layer makes for faces.
#
# `showlabel` survives, with a different meaning. Under Plots it set the `label` of the first two
# series, and the vector recipe had to pass `showlabel = i == 1` so that n paths did not put n copies
# of `"e-"` in the legend. A Makie legend is built by the wrapper from explicit elements (plan §1.3),
# so there is exactly one entry per species by construction and `showlabel` now only decides whether
# the legend is drawn at all.

"""
    driftpathplot(paths; kwargs...)
    driftpathplot!(ax, paths; kwargs...)

Plot the trajectories of the charge carriers of one `EHDriftPath`, or of a vector of them as returned
by `Event.drift_paths`.

The non-mutating form returns a `FigureAxisPlot` carrying the axis labels and the `e-` / `h+` legend;
the mutating form draws into an existing axis and sets neither. `plot(paths)` routes here.

Coordinates are drawn in millimetres, so a drift path can be drawn into the same axis as a
[`detectorplot`](@ref) of the detector it was drifted in.

## Keywords
* `showlabel`: draw the legend. Set it to `false` when overlaying paths on a figure that already has
  one.
* `electron_color`, `hole_color`: the colour of each species.
* `linewidth`, `linestyle`: the trajectories.
* `markersize`: size of the marker at each carrier's final position.
"""
@recipe DriftPathPlot (paths,) begin
    "Draw the `e-` / `h+` legend; acted on by `driftpathplot`, which owns the legend."
    showlabel = true
    "Colour of the electron trajectories."
    electron_color = Makie.automatic
    "Colour of the hole trajectories."
    hole_color = Makie.automatic
    "Width of the trajectories."
    linewidth = Makie.automatic
    "Style of the trajectories."
    linestyle = :solid
    "Size of the marker at each carrier's final position."
    markersize = Makie.automatic
end

const DriftPathLike = Union{EHDriftPath, AbstractVector{<:EHDriftPath}}

for GT in (:EHDriftPath, :(AbstractVector{<:EHDriftPath}))
    @eval Makie.plottype(::$GT) = DriftPathPlot
    @eval Makie.convert_arguments(::Type{<:DriftPathPlot}, d::$GT) = (d,)
end

Makie.preferred_axis_type(::DriftPathPlot) = Makie.Axis3


## ----------------------------------------------------------------------------------------------
## The recipe
## ----------------------------------------------------------------------------------------------

_drift_paths(dp::EHDriftPath) = [dp]
_drift_paths(dps::AbstractVector{<:EHDriftPath}) = dps

# The trajectories of one species, end to end, separated by a non-finite vertex so that Makie lifts
# the pen between paths instead of drawing a line from where one carrier was collected to where the
# next one started.
function _joined_paths(paths, species)
    pts = Makie.Point3f[]
    for dp in paths
        segment = point_coordinates(getfield(dp, species))
        isempty(segment) && continue
        isempty(pts) || push!(pts, Makie.Point3f(NaN32, NaN32, NaN32))
        append!(pts, segment)
    end
    pts
end

_final_positions(paths, species) =
    Makie.Point3f[point_coordinates(last(getfield(dp, species)))[1]
                  for dp in paths if !isempty(getfield(dp, species))]

function Makie.plot!(p::DriftPathPlot)
    paths = _drift_paths(p[1][])
    lw = p.linewidth[] === Makie.automatic ? SSD_LINEWIDTHS.driftpath : p.linewidth[]
    ms = p.markersize[] === Makie.automatic ? SSD_MARKERSIZES.drift_endpoint : p.markersize[]
    for (species, color) in ((:e_path, color_or_default(p.electron_color[], SSD_COLORS.electron)),
                             (:h_path, color_or_default(p.hole_color[], SSD_COLORS.hole)))
        line = _joined_paths(paths, species)
        isempty(line) && continue
        Makie.lines!(p, line; color, linewidth = lw, linestyle = p.linestyle[])
        Makie.scatter!(p, _final_positions(paths, species);
                       color, markersize = ms, strokewidth = 0)
    end
    pad_flat_extent!(p)
    return p
end


## ----------------------------------------------------------------------------------------------
## Figure-level wrapper
## ----------------------------------------------------------------------------------------------

function driftpathplot(gp::Makie.GridPosition, paths::DriftPathLike;
                       axis = (;), legend = (;), kwargs...)
    # A drift path is a curve through the detector, so it wants the axis a piece of detector geometry
    # wants: millimetres on all three dimensions and equal data scaling.
    gl = panel_layout(gp)
    ax = ssd_axis(gl[1, 1]; geometry_axis_setup().axis..., axis...)
    plt = driftpathplot!(ax, paths; kwargs...)
    finish_geometry_axis!(ax, plt, axis)
    if get(kwargs, :showlabel, true) && legend !== false
        add_legend!(gl,
                    drift_legend_entries(
                        electron_color = color_or_default(get(kwargs, :electron_color, Makie.automatic),
                                                        SSD_COLORS.electron),
                        hole_color = color_or_default(get(kwargs, :hole_color, Makie.automatic),
                                                    SSD_COLORS.hole))...;
                    corner = legend_corner(ax, plt), legend...)
    end
    Makie.AxisPlot(ax, plt)
end

function driftpathplot(paths::DriftPathLike; figure = (;), kwargs...)
    fig = Makie.Figure(; figure...)
    ap = driftpathplot(fig[1, 1], paths; kwargs...)
    Makie.FigureAxisPlot(fig, ap.axis, ap.plot)
end

Makie.plot(paths::DriftPathLike; kwargs...) = driftpathplot(paths; kwargs...)
