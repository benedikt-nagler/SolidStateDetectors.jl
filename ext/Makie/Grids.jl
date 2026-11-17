# Makie recipe and figure-level wrapper for the grid-point density of a `Grid`.
# 2 of the ~142 parity cases, and the first that needs a `GridLayout`.
#
# The Plots recipe asks for `layout --> length(iaxs)` and emits one `:stephist` series per grid axis,
# and Plots builds the subplots for it. A Makie recipe cannot create an axis, let alone several, so
# the split lands where plan D1 puts it: `GridPlot` draws the histogram of *one* axis into whatever
# axis it is given, and `gridplot` owns the layout, the panel titles and the labels.
#
# The labels are where this file departs from the golden spec on purpose. The Plots recipe writes
# "Grid point density - Axis 1" into the *x guide*, over ticks that are really a coordinate in
# internal metres or radians, and leaves the title and the y guide empty. Here that string becomes
# the panel title, x is labelled with the coordinate in display units like every other plot, and y
# names the count. Plan A1 ("readable units") is the rule that licenses it; `test/parity/gate.jl`
# carries the bridge, and the plan records it as the grids deviation.

# The count a bin of the histogram holds.
const GRID_COUNT_NAME = "Grid points"

"""
    gridplot(grid; kwargs...)
    gridplot!(ax, grid; dim, kwargs...)

Plot how densely the ticks of `grid` are spaced along each of its axes.

The non-mutating form returns a `Figure` with one panel per axis that has more than one tick, not a
`FigureAxisPlot`: a `GridLayout` has no single axis to hand back. The mutating form draws the
histogram of one axis, selected by `dim`, into an existing axis.

Use it to check where the refinement put its points before reading a potential off the same grid.

## Keywords
* `dim`: index into `grid.axes` of the axis to show. Mutating form only.
* `bins`: number of histogram bins. `automatic` is half the number of ticks, as in the Plots recipe.
* `color`, `linewidth`: the step outline.
"""
@recipe GridPlot (grid,) begin
    "Index into `grid.axes` of the axis to show."
    dim = 1
    "Number of histogram bins; `automatic` is half the number of ticks."
    bins = Makie.automatic
    "Colour of the step outline."
    color = Makie.automatic
    "Width of the step outline."
    linewidth = Makie.automatic
end

Makie.plottype(::Grid) = GridPlot
Makie.convert_arguments(::Type{<:GridPlot}, g::Grid) = (g,)

# Names and units of the three axes of a grid, by coordinate system. The angular axis is the one
# reason this cannot be a single method: it is the only grid axis not measured in metres.
grid_axis_name(::Grid{<:Any, 3, Cylindrical}, dim::Integer) = ("r", "φ", "z")[dim]
grid_axis_name(::Grid{<:Any, 3, Cartesian}, dim::Integer) = ("x", "y", "z")[dim]
grid_axis_unit(::Grid{<:Any, 3, Cylindrical}, dim::Integer) =
    dim == 2 ? internal_angle_unit : internal_length_unit
grid_axis_unit(::Grid{<:Any, 3, Cartesian}, ::Integer) = internal_length_unit

# The axes worth a panel. An axis with a single tick is a two-dimensional simulation's degenerate
# third dimension, and a one-bin histogram of it says nothing.
grid_plot_dims(grid::Grid) = findall(ax -> length(ax) > 1, collect(grid.axes))


## ----------------------------------------------------------------------------------------------
## The recipe
## ----------------------------------------------------------------------------------------------

"""
    grid_point_histogram(grid, dim, bins) -> (edges, counts)

Bin the ticks of axis `dim` of `grid`, with the edges rescaled to display units.

`bins` counts bins; `automatic` is `div(nticks, 2)`, which is what the Plots recipe asks for.
`histrange` then turns that into the nice bin width Plots also lands on, so the two backends put a
tick in the same bin.

The counts are taken on the ticks as the grid holds them, and only the edges are rescaled. Binning
the rescaled ticks instead moves a tick sitting exactly on a bin edge across it, and the two
backends then disagree by a count.
"""
function grid_point_histogram(grid::Grid, dim::Integer, bins = Makie.automatic)
    ticks = collect(grid[dim].ticks)
    nbins = bins === Makie.automatic ? max(div(length(ticks), 2), 1) : Int(bins)
    edges = StatsBase.histrange(ticks, nbins, :left)
    h = StatsBase.fit(StatsBase.Histogram, ticks, edges; closed = :left)
    (display_values(collect(edges), grid_axis_unit(grid, dim)), h.weights)
end

# The outline of a histogram as one polyline: up from zero at the first edge, along the top of every
# bin, and back down to zero at the last. Drawing it by hand rather than with `stephist!` is what
# keeps the bin edges the ones `grid_point_histogram` chose.
function step_outline(edges, counts)
    pts = [Makie.Point2f(edges[1], 0)]
    for (i, c) in enumerate(counts)
        push!(pts, Makie.Point2f(edges[i], c), Makie.Point2f(edges[i + 1], c))
    end
    push!(pts, Makie.Point2f(edges[end], 0))
    pts
end

function Makie.plot!(p::GridPlot)
    edges, counts = grid_point_histogram(p[1][], p.dim[], p.bins[])
    Makie.lines!(p, step_outline(edges, counts);
                 color = color_or_default(p.color[], SSD_PALETTE[1]),
                 linewidth = color_or_default(p.linewidth[], SSD_LINEWIDTHS.waveform))
    return p
end


"""
    count_ticks(maxcount, target = 5) -> Vector{Int}

Ticks at whole numbers from 0 to `maxcount`, roughly `target` of them.

The step is the first of 1, 2, 5, 10, ... that keeps the count of ticks at or below `target`. Makie's
automatic ticks are chosen for a continuous axis and put `2.5` on one whose values are counts, which
reads as a fraction of a grid point.
"""
function count_ticks(maxcount::Integer, target::Int = 5)
    maxcount <= 0 && return [0]
    for step in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
        maxcount / step <= target && return collect(0:step:maxcount)
    end
    collect(range(0, Int(maxcount), length = target))
end


## ----------------------------------------------------------------------------------------------
## Figure-level wrapper
## ----------------------------------------------------------------------------------------------

# Where panel `i` of `n` goes. Up to two panels share a row. A third wraps onto a second row, since
# three abreast leave each too narrow for its title.
grid_panel_cell(i::Int, n::Int) = n <= 2 ? (1, i) : (cld(i, 2), 2 - isodd(i))

function gridplot(gp::Makie.GridPosition, grid::Grid; axis = (;), kwargs...)
    gl = panel_layout(gp)
    dims = grid_plot_dims(grid)
    isempty(dims) && (dims = [1])
    axes_made = Makie.Axis[]
    for (col, dim) in enumerate(dims)
        _, counts = grid_point_histogram(grid, dim, get(kwargs, :bins, Makie.automatic))
        ax = ssd_axis(gl[grid_panel_cell(col, length(dims))...]; title = grid_panel_title(dim),
                      xname = grid_axis_name(grid, dim), yname = GRID_COUNT_NAME,
                      xunit = grid_axis_unit(grid, dim),
                      yticks = count_ticks(maximum(counts; init = 0)), axis...)
        gridplot!(ax, grid; dim, kwargs...)
        push!(axes_made, ax)
    end
    gl
end

function gridplot(grid::Grid; figure = (;), axis = (;), kwargs...)
    dims = grid_plot_dims(grid)
    isempty(dims) && (dims = [1])
    fig = Makie.Figure(; figure...)
    for (col, dim) in enumerate(dims)
        # The histogram is rebinned here only to reach its tallest bin, which is what fixes the tick
        # positions. Binning a few dozen ticks costs nothing next to drawing the panel.
        _, counts = grid_point_histogram(grid, dim, get(kwargs, :bins, Makie.automatic))
        ax = ssd_axis(fig[grid_panel_cell(col, length(dims))...]; title = grid_panel_title(dim),
                      xname = grid_axis_name(grid, dim), yname = GRID_COUNT_NAME,
                      xunit = grid_axis_unit(grid, dim),
                      yticks = count_ticks(maximum(counts; init = 0)), axis...)
        gridplot!(ax, grid; dim, kwargs...)
    end
    fig
end

Makie.plot(grid::Grid; kwargs...) = gridplot(grid; kwargs...)
