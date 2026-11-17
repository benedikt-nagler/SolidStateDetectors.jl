# Legend construction for the figure-level wrappers (plan §1.3).
#
# Plots builds the legend out of every series' `label`, and this package exploits that with a trick:
# each detector part emits an **empty `:path` series** purely so that the legend gets an entry drawn
# at full opacity while the geometry itself is drawn at `fillalpha = 0.2`.
#
# Makie collects legend entries only from plots that are *direct children of the axis*, so a
# `detectorplot` recipe that draws all contacts inside itself produces exactly zero entries. The
# figure-level wrapper therefore has to build the `Legend` from explicit elements — which is also
# strictly better: the swatch opacity is decoupled from the geometry opacity by construction, so the
# empty-series hack disappears rather than being ported.
#
# The label *strings* are not defined here. They come from `plot_label` in
# `src/PlotRecipes/plot_utils.jl`, which both backends share, so a Plots legend and a Makie legend
# can never drift apart.

"""
    legend_entries(det::SolidStateDetector; show_semiconductor = false, show_passives = true)

The legend elements and labels for a detector, in the order the Plots recipe emits them: contacts
first, then the semiconductor if shown, then the passives.

Returns `(elements, labels)` ready for `Makie.Legend(fig[…], elements, labels)`.
"""
function legend_entries(det::SolidStateDetector; seriestype::Symbol = :csg, kwargs...)
    # `detector_parts` (in `SolidStateDetector.jl`) is what fixes the order, and it is the same list
    # the recipe draws from, so a legend can never name a part the picture does not show.
    parts = detector_parts(det; kwargs...)
    (Makie.LegendElement[representation_element(seriestype, part_color(p)) for p in parts],
     String[plot_label(p) for p in parts])
end

"""
    representation_element(seriestype, color) -> LegendElement

Give the legend swatch matching how the geometry is drawn.

Mesh-based representations get a filled patch, `:wireframe` a line, `:samplesurface` a marker, and
`:slice`, whose samples trace an outline, a line.
"""
function representation_element(seriestype::Symbol, color)
    seriestype === :wireframe &&
        return Makie.LineElement(color = color, linewidth = SSD_LINEWIDTHS.wireframe)
    # A slice is a cloud of samples dense enough to read as an outline, so it is listed as one.
    seriestype === :slice && return Makie.LineElement(color = color, linewidth = 2)
    seriestype === :samplesurface &&
        return Makie.MarkerElement(color = color, marker = :circle, markersize = 10)
    geometry_legend_element(color)
end

"""
    geometry_legend_element(color) -> PolyElement

A legend swatch for a piece of geometry: opaque, with the dark outline the meshes are drawn with, no
matter how translucent the geometry itself is.
"""
geometry_legend_element(color) =
    Makie.PolyElement(color = color, strokecolor = :black, strokewidth = 1)

"""
    drift_legend_entries(; electron_color, hole_color, linewidth) -> (elements, labels)

The `"e-"` / `"h+"` legend of a drift-path plot, matching the Plots recipe's two labelled series.

The colours are arguments rather than constants so that the legend cannot describe one pair of
colours while `driftpathplot` draws another.
"""
drift_legend_entries(; electron_color = SSD_COLORS.electron, hole_color = SSD_COLORS.hole,
                     linewidth::Real = SSD_LINEWIDTHS.driftpath) = (
    Makie.LegendElement[Makie.LineElement(color = electron_color, linewidth = linewidth),
                        Makie.LineElement(color = hole_color, linewidth = linewidth)],
    String[electron_label, hole_label],
)

"""
    eh_contribution_legend_entries(colors; linewidth) -> (elements, labels)

The three-entry legend of an electron/hole contribution plot.

`colors` is the named tuple `_eh_colors` resolved, so the legend cannot describe one set of colours
while `ehcontributionplot` draws another.
"""
eh_contribution_legend_entries(colors; linewidth::Real = SSD_LINEWIDTHS.contribution) = (
    Makie.LegendElement[Makie.LineElement(color = colors.electron, linewidth = linewidth),
                        Makie.LineElement(color = colors.hole, linewidth = linewidth),
                        Makie.LineElement(color = colors.sum, linewidth = linewidth)],
    String[electron_contribution_label, hole_contribution_label, sum_contribution_label],
)

"""
    waveform_legend_entries(labels; color, linewidth) -> (elements, labels)

One legend entry per waveform, in the colours `waveformplot` draws them in.

`labels` names the waveforms in order, normally by contact. An unlabelled waveform plot gets no
legend at all: Plots fills one with its auto labels `y1`, `y2`, which name nothing.
"""
function waveform_legend_entries(labels; color = Makie.automatic,
                                 linewidth::Real = SSD_LINEWIDTHS.waveform)
    strs = String[string(l) for l in labels]
    (Makie.LegendElement[Makie.LineElement(color = _waveform_color(color, i), linewidth = linewidth)
                         for i in eachindex(strs)],
     strs)
end

"""
    charge_cloud_legend_element(color) -> MarkerElement

A legend swatch for a charge cloud: one carrier, drawn opaque at a readable size whatever the cloud
itself is drawn at.
"""
charge_cloud_legend_element(color) =
    Makie.MarkerElement(color = color, marker = :circle, markersize = 10)

"""
    add_legend!(figure, elements, labels; position = (1, 2), kwargs...)

Place a `Legend` built from explicit elements, doing nothing when there is nothing to show.
"""
function add_legend!(target, elements, labels;
                     position::Tuple{Int, Int} = (1, 1), corner::Symbol = :rt, kwargs...)
    isempty(elements) && return nothing
    # No corner free of data: the legend goes beside the axis instead of on top of it.
    corner === :outside && (position = (1, 2))
    # Any cell but the axis' own is a legend beside the plot, sized by the layout as usual.
    position == (1, 1) || return Makie.Legend(target[position...], elements, labels; kwargs...)
    halign, valign = LEGEND_CORNERS[corner]
    Makie.Legend(target[1, 1], elements, labels;
                 tellwidth = false, tellheight = false, halign, valign,
                 margin = (8, 8, 8, 8), padding = (6, 6, 4, 4), patchsize = (18, 12),
                 rowgap = 1, backgroundcolor = (:white, 0.85), kwargs...)
end

const LEGEND_CORNERS = (rt = (:right, :top), lt = (:left, :top),
                        rb = (:right, :bottom), lb = (:left, :bottom))

"""
    best_legend_corner(plot) -> Symbol

The corner of a two-dimensional axis a legend fits into without covering data: `:rt`, `:lt`, `:rb`
or `:lb`, tried in that order, or `:outside` when every corner has data in it.
"""
function best_legend_corner(plot; width::Real = 0.42, height::Real = 0.3)
    pts = Makie.Point2f[]
    _collect_points!(pts, plot)
    filter!(pt -> all(isfinite, pt), pts)
    isempty(pts) && return :rt
    xs, ys = first.(pts), last.(pts)
    (x0, x1), (y0, y1) = extrema(xs), extrema(ys)
    (x1 > x0 && y1 > y0) || return :rt
    u = (xs .- x0) ./ (x1 - x0)
    v = (ys .- y0) ./ (y1 - y0)
    hits(right, top) = count(i -> (right ? u[i] > 1 - width : u[i] < width) &&
                                  (top ? v[i] > 1 - height : v[i] < height), eachindex(u))
    corners = (:rt => hits(true, true), :lt => hits(false, true),
               :rb => hits(true, false), :lb => hits(false, false))
    i = findfirst(c -> last(c) == 0, corners)
    i === nothing ? :outside : first(corners[i])
end

function _collect_points!(pts, plot)
    children = plot.plots
    if isempty(children)
        try
            for pt in Makie.point_iterator(plot)
                push!(pts, Makie.Point2f(pt[1], pt[2]))
            end
        catch
        end
    else
        foreach(c -> _collect_points!(pts, c), children)
    end
    pts
end

# Where a legend goes inside its axis. A three-dimensional axis is fitted to its bounding sphere,
# which leaves its top corners free, and the z ticks run down the right side. A flat axis is searched
# for a corner without data.
legend_corner(::Makie.Axis3, plot) = :lt
legend_corner(::Any, plot) = best_legend_corner(plot)
