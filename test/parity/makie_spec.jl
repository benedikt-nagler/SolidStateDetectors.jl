# Tier-1 of the parity harness: turn a Makie figure into the same backend-neutral `CaseSpec` that
# `devscripts/plots_spec.jl` extracts from Plots.jl, so the two can be diffed field by field.
#
# What the two sides deliberately do *not* share is where information lives:
#
#   * Plots keeps title, guides, limits, `clims`, colorbar and legend on the subplot; a Makie recipe
#     cannot set any of them (plan §1.1), so on this side they come from the `Axis`, the `Colorbar`
#     and the `Legend` blocks of the figure — i.e. from the figure-level wrapper functions. A recipe
#     tested on its own will show them empty, which is the point.
#   * A Plots heatmap `z` is `(length(y), length(x))`; a Makie one is `(length(x), length(y))`
#     (plan §1.7). Both are recorded as the backend holds them; `compare_data(...; transposed = true)`
#     is what bridges them.
#
# Requires a Makie backend to be loaded (CairoMakie in the test environment).

module MakieSpec

using Makie

include(joinpath(@__DIR__, "spec.jl"))
using .PlotParitySpec

export makie_case_spec

const AXIS_TYPES = Union{Makie.Axis, Makie.Axis3, Makie.PolarAxis}

# Plot types whose `colormap` attribute is actually in use, i.e. that map numbers to colours. A
# `mesh!` or `scatter!` drawn in one flat colour carries the attribute too and must not report it.
const COLORMAPPED_PLOTKEYS = (:heatmap, :image, :surface, :volume, :contourf, :contour)

## ----------------------------------------------------------------------------------------------
## Attribute readers
## ----------------------------------------------------------------------------------------------

# Blocks differ in which attributes they have at all: a PolarAxis has no `xlabel`, an Axis no
# `zlabel`. Missing means "" / nothing rather than an error.
function _attr(obj, key, default = nothing)
    v = try
        Makie.to_value(getproperty(obj, key))
    catch
        return default
    end
    v === nothing ? default : v
end

_str(::Nothing) = ""
_str(x) = string(x)
_str(x::Symbol) = String(x)

_num(x::Real) = Float64(x)
_num(::Any) = nothing

_hex2(x) = string(round(Int, clamp(x, 0, 1) * 255), base = 16, pad = 2)
function _hex(c)
    rgba = Makie.to_color(c)
    "#" * _hex2(rgba.r) * _hex2(rgba.g) * _hex2(rgba.b) * (rgba.alpha < 1 ? _hex2(rgba.alpha) : "")
end

# Colormaps are recorded as three sampled stops, exactly as on the Plots side, so `:inferno` on one
# side and a `ContinuousColorGradient` on the other produce comparable strings. Use
# `PlotParitySpec.compare_colormap` to diff them — the two backends interpolate slightly differently.
function _colormap_string(cmap)
    cmap === nothing && return ""
    try
        cs = Makie.to_colormap(cmap)
        isempty(cs) && return ""
        stops = [cs[1], cs[cld(length(cs), 2)], cs[end]]
        return "grad(" * join(_hex.(stops), ",") * ")"
    catch
        return _str(cmap)
    end
end

# A single colour, a colour per vertex, or a numeric array mapped through the colormap.
function _color_string(c)
    c === nothing && return ""
    c isa Symbol && return String(c)
    c isa AbstractArray && return "array($(length(c)))"
    try
        return _hex(c)
    catch
        return _str(c)
    end
end

_colorrange(v) = (v === nothing || v isa Makie.Automatic) ? nothing :
                 (v isa Tuple || v isa AbstractVector) && length(v) == 2 ?
                     Float64[Float64(v[1]), Float64(v[2])] : nothing

# The golden spec records the limits Plots was *asked* for (`sp[:xaxis][:lims]`), not the ones GR
# renders, so the comparable quantity here is `ax.limits[]` — `finallimits` additionally carries the
# expansion `DataAspect()` applies to fill the axis box, which GR does at render time too. A
# dimension the caller left open falls back to `finallimits`, which is the only value there is.
function _limits(ax)
    fallback = let r = _attr(ax, :finallimits)
        r === nothing ? (nothing, nothing) :
            (Float64[r.origin[1], r.origin[1] + r.widths[1]],
             Float64[r.origin[2], r.origin[2] + r.widths[2]])
    end
    req = _attr(ax, :limits)
    req === nothing && return fallback
    dims = req isa Tuple && length(req) == 4 ? ((req[1], req[2]), (req[3], req[4])) :
           req isa Tuple && length(req) == 2 ? (req[1], req[2]) : nothing
    dims === nothing && return fallback
    ntuple(2) do i
        d = dims[i]
        (d isa Tuple && length(d) == 2 && all(x -> x isa Real, d)) ?
            Float64[d[1], d[2]] : fallback[i]
    end
end

_aspect_string(a) = a === nothing ? "" : a isa Makie.DataAspect ? "1" : _str(a)

## ----------------------------------------------------------------------------------------------
## Walking the plot tree
## ----------------------------------------------------------------------------------------------

# Makie's atomic plot types: the ones that carry geometry themselves. Everything else is a recipe,
# which owns only child plots.
const ATOMIC_PLOTKEYS = (:lines, :linesegments, :scatter, :meshscatter, :mesh, :surface,
                         :heatmap, :image, :volume, :voxels, :text)

# The two invisible points `pad_flat_extent!` adds so that a flat object cannot collapse an axis.
# They widen the limits and draw nothing, so they are not series data.
function is_extent_anchor(p)
    Makie.plotkey(p) === :scatter || return false
    ms = _attr(p, :markersize, nothing)
    ms !== nothing && all(iszero, ms)
end

"""
    leaf_plots(p) -> Vector

The drawable primitives under `p`. A recipe plot is a *composite*: it owns no geometry itself, only
child plots, so a parity comparison against Plots' flat `series_list` has to compare against the
leaves.

A recipe that drew **nothing** has no children either, and is not a primitive: it contributes no
series. Returning it would read its first positional argument as series data, which for
`fieldlines_car_x` is a whole `Simulation`.
"""
function leaf_plots(p)
    children = _attr(p, :plots, nothing)
    if children === nothing || isempty(children)
        return (Makie.plotkey(p) in ATOMIC_PLOTKEYS && !is_extent_anchor(p)) ? [p] : []
    end
    reduce(vcat, (leaf_plots(c) for c in children); init = [])
end

# `mesh!` carries its faces in the geometry rather than in an attribute.
function _n_connections(p)
    Makie.plotkey(p) === :mesh || return 0
    try
        return length(Makie.GeometryBasics.faces(p[1][]))
    catch
        return 0
    end
end

"""
    polar_cell_matrix(p) -> Matrix or nothing

The values of a polar cross section, as an (angular, radial) matrix, read back off the cell mesh
`polar_cells!` draws; `nothing` for any other plot.

The mesh holds one quad per cell, four vertices of one value each, with the radial index running
fastest. An angular cell wider than `POLAR_MAX_CELL` is split into equal parts and so appears as
several identical rows.
"""
function polar_cell_matrix(p)
    Makie.plotkey(p) === :mesh || return nothing
    c = _attr(p, :color)
    (c isa AbstractVector{<:Real} && length(c) >= 4 && length(c) % 4 == 0) || return nothing
    pts = Makie.GeometryBasics.coordinates(p[1][])
    n_cells = length(c) ÷ 4
    θ0 = pts[1][1]
    n_radial = something(findfirst(k -> pts[4k + 1][1] != θ0, 1:(n_cells - 1)), n_cells)
    permutedims(reshape(c[1:4:end], n_radial, :))
end

# Positional arguments of a primitive, as the backend converted them.
function _arg(p, i)
    try
        return p[i][]
    catch
        return nothing
    end
end

"""
    _xyz(p) -> (x, y, z)

The positional data of a primitive, in the shape the Plots golden spec holds it.

Point-based primitives (`lines`, `scatter`, `linesegments`, and the `lines` a `contour!` decomposes
into) carry a single `Vector{Point{N}}` where Plots carried separate `x`, `y` and `z` vectors, so
that one argument is split back into its components. Everything else — a heatmap's `(x, y, z)`, for
instance — is passed through.
"""
function _xyz(p)
    a = _arg(p, 1)
    if a isa AbstractVector && !isempty(a) && eltype(a) <: AbstractVector{<:Real}
        n = length(first(a))
        return ntuple(i -> i <= n ? [q[i] for q in a] : nothing, 3)
    end
    (a, _arg(p, 2), _arg(p, 3))
end

function _series_spec(p, subplot_idx::Int)
    x, y, z = _xyz(p)
    # A `surface!` used as a polar pcolor carries zero heights in its third argument and the numbers
    # in `color` (the idiom Makie's own PolarAxis documentation uses), so the colour-mapped matrix is
    # what has to be recorded as `z` — that is what the Plots heatmap `z` is compared against.
    c = _attr(p, :color)
    c isa AbstractMatrix{<:Real} && (z = c)
    # A polar cross section is a mesh of flat cells: its values are recorded as the field `z`, and
    # the positions, which are cell corners, are not series data.
    cells = polar_cell_matrix(p)
    cells === nothing || ((x, y, z) = (nothing, nothing, cells))
    SeriesSpec(
        subplot_idx,
        cells === nothing ? _str(Makie.plotkey(p)) : "polarcells",
        _str(_attr(p, :label, "")),
        _color_string(_attr(p, :color)),
        _color_string(_attr(p, :strokecolor, _attr(p, :color))),
        _num(_attr(p, :linewidth, _attr(p, :strokewidth))),
        _num(_attr(p, :alpha)),
        _num(_attr(p, :alpha)),
        _str(_attr(p, :marker, "")),
        _num(_attr(p, :markersize)),
        _n_connections(p),
        data_spec(x), data_spec(y), data_spec(z),
    )
end

# The axis unit is not an axis attribute in Makie the way it is in Plots: with a
# `UnitfulConversion` installed the units ride on the data. Reading it back off the first leaf that
# has data for that dimension is therefore the faithful equivalent of Plots' `sp[:xaxis][:unit]`.
function _axis_unit(series, dim::Symbol)
    for s in series
        d = getfield(s, dim)
        d === nothing || isempty(d.unit) || return d.unit
    end
    ""
end

function _axis_spec(ax, series, colorbar, legend)
    xlims, ylims = _limits(ax)
    cr = nothing
    for p in series_plots_of(ax)
        cr = _colorrange(_attr(p, :colorrange))
        cr === nothing || break
    end
    # Only from a plot that actually maps values through the colormap. Every Makie plot type carries
    # a `colormap` attribute (viridis by default), so scanning for a non-`nothing` one would give a
    # geometry plot a colormap it never uses, where the Plots side records the colormap of a
    # *gradient-coloured series* and leaves it empty for a mesh.
    cmap = ""
    for p in series_plots_of(ax)
        (Makie.plotkey(p) in COLORMAPPED_PLOTKEYS || polar_cell_matrix(p) !== nothing) || continue
        c = _attr(p, :colormap)
        c === nothing && continue
        cmap = _colormap_string(c)
        break
    end
    AxisSpec(
        _str(_attr(ax, :title, "")),
        _str(_attr(ax, :xlabel, "")), _str(_attr(ax, :ylabel, "")), _str(_attr(ax, :zlabel, "")),
        _axis_unit(series, :x), _axis_unit(series, :y), _axis_unit(series, :z),
        xlims, ylims, cr,
        colorbar === nothing ? "none" : "colorbar",
        colorbar === nothing ? "" : _str(_attr(colorbar, :label, "")),
        cmap,
        _aspect_string(_attr(ax, :aspect)),
        ax isa Makie.PolarAxis ? "polar" : ax isa Makie.Axis3 ? "3d" : "none",
        legend === nothing ? "none" : _str(_attr(legend, :position, "best")),
        legend === nothing ? String[] : _legend_entries(legend),
    )
end

series_plots_of(ax) = reduce(vcat, (leaf_plots(p) for p in ax.scene.plots); init = [])

function _legend_entries(legend)
    entries = String[]
    try
        for group in legend.entrygroups[]
            for entry in group[2]
                push!(entries, string(entry.label[]))
            end
        end
    catch
    end
    entries
end

## ----------------------------------------------------------------------------------------------
## Figure -> CaseSpec
## ----------------------------------------------------------------------------------------------

_figure_of(x::Makie.Figure) = x
_figure_of(x::Makie.FigureAxisPlot) = x.figure
_figure_of(x) = error("cannot extract a Figure from a $(typeof(x))")

"""
    makie_case_spec(figlike; id, group, call, note = "") -> CaseSpec

Introspect a `Figure` or `FigureAxisPlot` into a spec comparable with the Plots golden spec.

`update_state_before_display!` is called first: until the figure has been laid out, `finallimits`
still holds the placeholder `(0, 10) × (0, 10)` and every limit assertion would compare against
nonsense.
"""
function makie_case_spec(figlike; id::AbstractString, group::AbstractString,
                         call::AbstractString, note::AbstractString = "")
    fig = _figure_of(figlike)
    Makie.update_state_before_display!(fig)

    axes_blocks = [b for b in fig.content if b isa AXIS_TYPES]
    colorbars = [b for b in fig.content if b isa Makie.Colorbar]
    legends = [b for b in fig.content if b isa Makie.Legend]

    series = SeriesSpec[]
    axspecs = AxisSpec[]
    for (i, ax) in enumerate(axes_blocks)
        own = SeriesSpec[_series_spec(p, i) for p in series_plots_of(ax)]
        append!(series, own)
        # With a single axis — the overwhelming majority of cases — any colorbar and legend in the
        # figure belong to it. With several, they are matched by order of creation.
        push!(axspecs, _axis_spec(ax, own,
                                  get(colorbars, i, length(colorbars) == 1 ? colorbars[1] : nothing),
                                  get(legends, i, length(legends) == 1 ? legends[1] : nothing)))
    end

    CaseSpec(; id, group, call, backend = "makie", axes = axspecs, series, note)
end

end # module
