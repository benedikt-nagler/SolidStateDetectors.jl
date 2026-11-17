# Axis construction for the figure-level wrappers.
#
# A Makie recipe cannot set the axis title, the labels, the aspect ratio or the limits (plan §1.1) —
# every one of those is a property of the `Axis`, which the *caller* creates. Since a Plots recipe
# sets all four routinely, reproducing them is the whole job of the figure-level wrapper functions,
# and `ssd_axis` is the one place that does it.
#
# It also owns the two things that are easy to get subtly wrong:
#
#   * the unit conversions (D2) — `xunit`/`yunit` are the units of the *data* the recipe will pass,
#     and the axis converts them to `display_unit` for the ticks while `unit_label` writes the unit
#     into the label;
#   * limits — with a conversion installed these live in *converted* space, so a limit given in
#     metres has to be pushed through `display_value` before it reaches the axis. Passing Unitful
#     limits and letting `ssd_axis` do the conversion is the safe path.

# `preferred_axis_type` runs before the recipe's default attributes are merged into the plot, so the
# only place the user-passed cross-section keywords exist yet is the raw kwargs dict `p.kw`. An absent
# key means "not selected" (missing) — exactly the default used by `get_crosssection_idx_and_value`.
_selector(p, key) = Makie.to_value(get(p.kw, key, missing))

"""
    ssd_axis(gridposition; kwargs...) -> Axis / PolarAxis / Axis3

Create the axis a figure-level wrapper draws into.

## Keywords
* `axistype`: `Makie.Axis` (default), `Makie.PolarAxis` for cylindrical `z` cross sections, or
  `Makie.Axis3` for 3-D geometry.
* `title`: axis title, the string the Plots recipes set via `title -->`.
* `xname`, `yname`, `zname`: quantity names (`"r"`, `"z"`, …). Combined with the units into
  `"r / mm"`-style labels. A `PolarAxis` has no labels and ignores them.
* `xunit`, `yunit`, `zunit`: the units of the data that will be plotted, in internal units
  (`internal_length_unit`, `internal_angle_unit`, …). `nothing` or `NoUnits` means no conversion.
* `xlims`, `ylims`, `zlims`: `(lo, hi)` tuples, either Unitful quantities or bare numbers already in
  the data's unit. Converted to display space here. `zlims` applies to an `Axis3` only.
* `aspect_ratio_one`: `true` reproduces the Plots recipes' `aspect_ratio --> 1` as `DataAspect()`.
* `dim_conversions`: install `UnitfulConversion`s so that **Unitful** data can be plotted into this
  axis. Off by default, because the recipes rescale their data with `display_values` instead: Makie
  0.24's `heatmap` cannot take Unitful data at all (it recurses until the stack overflows), and an
  axis carrying a conversion then refuses the bare numbers the heatmap does accept. Turn it on for
  an axis that will only receive point-based plots (`lines`, `scatter`) of Unitful data.

Any further keywords are forwarded to the axis constructor.
"""
function ssd_axis(gridposition;
                  axistype = Makie.Axis,
                  title::AbstractString = "",
                  xname::AbstractString = "", yname::AbstractString = "", zname::AbstractString = "",
                  xunit = nothing, yunit = nothing, zunit = nothing,
                  xlims = nothing, ylims = nothing, zlims = nothing,
                  aspect_ratio_one::Bool = false,
                  dim_conversions::Bool = false,
                  kwargs...)
    attrs = Dict{Symbol, Any}(:title => title)
    if axistype === Makie.PolarAxis
        # The title of an `Axis` is bold and that of a `PolarAxis` is not; one style for both, and
        # the same faint grid.
        attrs[:titlefont] = :bold
        attrs[:titlesize] = Makie.to_value(Makie.current_default_theme()[:fontsize])
        attrs[:spinewidth] = 1
        attrs[:rgridcolor] = (:black, 0.12)
        attrs[:thetagridcolor] = (:black, 0.12)
    end

    # A PolarAxis has neither labels nor a `limits` in the Cartesian sense, and forcing an aspect on
    # it is meaningless — it is round by construction.
    if axistype !== Makie.PolarAxis
        # The label carries the *display* unit — the axis shows millimetres, so the label must say
        # "r / mm" even though the data arrives in metres.
        attrs[:xlabel] = unit_label(xname, display_unit(xunit))
        attrs[:ylabel] = unit_label(yname, display_unit(yunit))
        if axistype === Makie.Axis3
            attrs[:zlabel] = unit_label(zname, display_unit(zunit))
            # Makie's default `:fitzoom` scales the axis to fit its *projected bounding box* and
            # then draws the decorations outside that, so an elongated object pushes the z label off
            # the figure entirely — silently, since the attribute is set correctly and only the
            # rendered image shows it missing. `:fit` fits the bounding sphere: a little more white
            # space, every label inside.
            attrs[:viewmode] = :fit
            # The camera of the Plots recipes (GR's `camera = (30, 30)`): x runs along the front
            # edge, y recedes to the right. Makie's default looks from the opposite quadrant, which
            # mirrors every detector picture in the manual.
            attrs[:azimuth] = -pi / 3
            attrs[:elevation] = pi / 6
            # Room for tick labels as long as "-1000" between the axis and its label.
            attrs[:zlabeloffset] = 58
            attrs[:ylabeloffset] = 44
            attrs[:xlabeloffset] = 40
            attrs[:protrusions] = 40
        end
        # `Axis3` spells equal data scaling as the symbol `:data`; `Axis` as `DataAspect()`.
        aspect_ratio_one && (attrs[:aspect] = axistype === Makie.Axis3 ? :data : Makie.DataAspect())
        lims = (_display_lims(xlims, xunit)..., _display_lims(ylims, yunit)...,
                (axistype === Makie.Axis3 ? _display_lims(zlims, zunit) : ())...)
        any(!isnothing, lims) && (attrs[:limits] = lims)
    end

    conv = dim_conversions ? _conversion_kwargs(axistype, xunit, yunit, zunit) : ()
    axistype(gridposition; attrs..., conv..., kwargs...)
end

_display_lims(::Nothing, _) = (nothing, nothing)
_display_lims(lims, u) = (_display_lim(lims[1], u), _display_lim(lims[2], u))
_display_lim(::Nothing, _) = nothing
_display_lim(v::Unitful.Quantity, _) = display_value(v)
_display_lim(v::Real, u) = u === nothing ? Float64(v) : display_value(v * u)

# `dim3_conversion` only exists on an Axis3; setting it elsewhere errors. Dimensions without a unit
# are left out entirely rather than set to a "no conversion" value.
function _conversion_kwargs(axistype, xunit, yunit, zunit)
    dims = axistype === Makie.Axis3 ?
        ((:dim1_conversion, xunit), (:dim2_conversion, yunit), (:dim3_conversion, zunit)) :
        ((:dim1_conversion, xunit), (:dim2_conversion, yunit))
    kw = Pair{Symbol, Any}[]
    for (key, u) in dims
        conv = unit_conversion(u)
        conv === nothing || push!(kw, key => conv)
    end
    kw
end

"""
    color_limits(data) -> (lo, hi)
    color_limits(lo, hi) -> (lo, hi)

The colour range Makie would resolve `Makie.automatic` to for `data`, as plain `Float64`s.

The figure-level wrappers have to work this out themselves: a `Colorbar` needs the range before the
figure is displayed while `plot.colorrange[]` is still `Automatic()` (plan §1.4). It therefore has to
agree with Makie's own resolution in every case, including the degenerate one — a *constant* slice.
There, an explicit `colorrange` with equal ends throws (`"Can't interpolate in a range where
cmin == cmax"`) while `automatic` on the same data renders, because Makie widens it by half a unit
(`distinct_extrema_nan`); the same widening is applied here. A constant slice is not exotic: the
effective charge density of an intrinsic detector is identically zero.
"""
function color_limits(lo::Real, hi::Real)
    l, h = Float64(lo), Float64(hi)
    l == h ? (l - 0.5, h + 0.5) : (l, h)
end
color_limits(range::Tuple{Real, Real}) = color_limits(range[1], range[2])
color_limits(data::AbstractArray) = color_limits(Makie.extrema_nan(data)...)

"""
    add_colorbar!(figure, plot; label = "", colorrange = nothing, position = (1, 2))

Attach the `Colorbar` a Plots recipe would have produced implicitly from `clims` + `colorbar_title`.

An explicit `colorrange` matters: until a figure is displayed, `plot.colorrange[]` is still
`Makie.automatic`, so a colorbar built from the plot alone can come out blank or wrongly scaled
(plan §1.4). `_get_potential_plot_information` already returns exactly the `clims` Plots used.

`plot` may be a composite recipe plot: the colorbar is built from its colour-mapped child rather
than from the composite itself, because a recipe that draws a heatmap *and* a contour overlay on top
of it carries two colormaps, and `Colorbar` refuses to guess between them.
"""
function add_colorbar!(target, plot;
                       label::AbstractString = "", colorrange = nothing, axis = nothing,
                       position::Union{Nothing, Tuple{Int, Int}} = nothing, kwargs...)
    attrs = Dict{Symbol, Any}(:label => label)
    colorrange === nothing || (attrs[:colorrange] = colorrange)
    source = colorbar_source(plot)
    # An explicit cell, or no axis to follow: an ordinary colorbar in the layout.
    (position === nothing && axis !== nothing) ||
        return Makie.Colorbar(target[something(position, (1, 2))...], source; attrs..., kwargs...)
    attached_colorbar!(target, axis, source; attrs..., kwargs...)
end

"""
    attached_colorbar!(layout, axis, source; gap = 10, width = 14, kwargs...) -> Colorbar

Place a colorbar against the right edge of what `axis` draws, as tall as the drawing.

An axis with equal data scaling is narrower or shorter than its layout cell. A colorbar in the
neighbouring cell then floats away from a narrow plot and towers over a flat one. This one follows
the plot area instead. Column 2 of `layout` is kept as wide as the colorbar and its labels need, so
nothing is drawn outside the figure.
"""
function attached_colorbar!(layout, axis, source; gap::Real = 10, width::Real = 14, kwargs...)
    bbox = Makie.lift(axis.scene.viewport) do vp
        x0, y0 = Makie.minimum(vp)
        w, h = Makie.widths(vp)
        if axis isa Makie.PolarAxis
            # The disc is inscribed in the viewport, with its angle labels just outside it.
            radius = min(w, h) / 2
            cx, cy = x0 + w / 2, y0 + h / 2
            Makie.Rect2f(cx + radius + 34, cy - 0.85 * radius, width, 1.7 * radius)
        else
            Makie.Rect2f(x0 + w + gap, y0, width, h)
        end
    end
    cb = Makie.Colorbar(axis.parent, source; bbox, kwargs...)
    reserve = Makie.lift(cb.layoutobservables.protrusions) do pr
        Makie.Fixed(gap + width + pr.right + (axis isa Makie.PolarAxis ? 24 : 0))
    end
    Makie.Box(layout[1, 2]; visible = false, width = Makie.lift(r -> r.x, reserve))
    cb
end

# The child whose colour scale the colorbar describes: the field itself, not an overlay drawn on top
# of it. Depth-first, so it also finds the heatmap of a recipe that wraps another recipe.
const COLORBAR_PLOTKEYS = (:heatmap, :image, :surface, :mesh, :volume, :contourf)

function colorbar_source(plot)
    children = Makie.to_value(getproperty(plot, :plots))
    (children === nothing || isempty(children)) && return plot
    for child in children
        Makie.plotkey(child) in COLORBAR_PLOTKEYS && return child
    end
    for child in children
        found = colorbar_source(child)
        found === child || return found
    end
    plot
end


"""
    panel_layout(gridposition) -> GridLayout

A nested layout for one plot placed into a caller's figure.

The axis goes in `[1, 1]` and its colorbar or legend in `[1, 2]`, so a plot dropped into one cell
keeps its own furniture instead of colliding with a neighbour's. This is what lets every figure-level
wrapper take a `GridPosition` as well as build its own `Figure`, and so what makes the multi-panel
figures of the manual possible.
"""
panel_layout(gp::Makie.GridPosition) = Makie.GridLayout(gp)


"""
    pad_flat_extent!(plot; minfrac = 0.2)

Keep a three-dimensional plot from collapsing an axis.

A single point, a flat surface or a thin plate has no extent in at least one dimension, and equal
data scaling then squeezes that dimension to nothing and stacks its tick labels on top of each other.
Two invisible anchor points widen every such dimension to `minfrac` of the largest extent (one
display unit for a lone point). Anchors live in the plot itself, so the axis keeps autoscaling when
more plots are added to it.
"""
function pad_flat_extent!(p; minfrac::Real = 0.2)
    bb = Makie.data_limits(p)
    lo, w = Float64.(collect(Makie.minimum(bb))), Float64.(collect(Makie.widths(bb)))
    (length(lo) == 3 && all(isfinite, lo) && all(isfinite, w)) || return p
    largest = maximum(w)
    need = largest > 0 ? minfrac * largest : 1.0
    any(<(need), w) || return p
    hi = lo .+ w
    for i in 1:3
        w[i] < need || continue
        c = lo[i] + w[i] / 2
        lo[i], hi[i] = c - need / 2, c + need / 2
    end
    Makie.scatter!(p, [Makie.Point3f(lo...), Makie.Point3f(hi...)];
                   markersize = 0, color = :transparent, strokewidth = 0, inspectable = false)
    p
end


"""
    finish_geometry_axis!(ax, plot, axis_kwargs)

Fit the ticks of a three-dimensional axis to the shape of what was drawn into it.

Under equal data scaling a flat object leaves one axis a fraction of the length of the others, and
five tick labels do not fit along it. Each dimension gets a number of ticks in proportion to its
extent, down to two. A dimension whose ticks the caller set through `axis_kwargs` is left alone, and
so is every axis that is not an `Axis3`.
"""
function finish_geometry_axis!(ax::Makie.Axis3, plt, axis_kwargs)
    bb = Makie.data_limits(plt)
    lo, w = Float64.(collect(Makie.minimum(bb))), Float64.(collect(Makie.widths(bb)))
    largest = maximum(w)
    (isfinite(largest) && largest > 0) || return ax
    for (i, key) in enumerate((:xticks, :yticks, :zticks))
        haskey(axis_kwargs, key) && continue
        n = clamp(round(Int, 5 * w[i] / largest), 2, 5)
        # A small object far from the origin has long tick labels: "1000.2", not "0.2".
        offset = maximum(abs, (lo[i], lo[i] + w[i]))
        w[i] > 0 && offset / w[i] > 50 && (n = min(n, 3))
        n < 5 && setproperty!(ax, key, Makie.LinearTicks(n))
    end
    ax
end
finish_geometry_axis!(ax, plt, axis_kwargs) = ax


"""
    readable_on(values, colormap, colorrange) -> Symbol

Give `:white` or `:black`, whichever reads better on top of `values` drawn with `colormap`.

The radial tick labels of a polar cross section sit on the data. White disappears on the bright end
of a colormap and black on the dark end, so the choice follows the mean luminance under the labels.
"""
function readable_on(values, colormap, colorrange)
    colors = Makie.to_colormap(colormap)
    lo, hi = colorrange
    finite = filter(isfinite, Float64.(vec(collect(values))))
    (isempty(finite) || hi <= lo) && return :white
    luminance = sum(finite) do v
        t = clamp((v - lo) / (hi - lo), 0, 1)
        c = colors[round(Int, 1 + t * (length(colors) - 1))]
        0.299 * Makie.red(c) + 0.587 * Makie.green(c) + 0.114 * Makie.blue(c)
    end / length(finite)
    luminance > 0.55 ? :black : :white
end
