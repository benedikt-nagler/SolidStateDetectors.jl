# Tier-0 of the parity harness: turn a Plots.jl `Plot` object into a backend-neutral `CaseSpec`.
#
# A Plots plot object can be fully introspected *without rendering* — no display, no GKS, no PNG —
# which is what makes the golden specs cheap to regenerate and independent of GR's output.
#
# The two structural facts this file encodes (both verified by probing, both load-bearing for the
# migration):
#
#   * Plots' Unitful integration strips units off the series data and stores them per axis
#     (`sp[:xaxis][:unit]`), leaving the guide bare. So the units go into the `AxisSpec`, and the
#     `DataSpec`s hold plain numbers.
#   * A heatmap/contour `z` is a `Plots.Surface` whose matrix is `(length(y), length(x))` — the
#     transpose of what Makie wants. It is stored here exactly as Plots holds it; undoing the
#     transpose is the *comparison's* job (`compare_data(...; transposed = true)`), so the golden
#     spec stays a faithful record of Plots' behaviour.

module PlotsSpec

using Plots
using Unitful

include(joinpath(@__DIR__, "..", "test", "parity", "spec.jl"))
using .PlotParitySpec

export plots_case_spec

## ----------------------------------------------------------------------------------------------
## Attribute readers
## ----------------------------------------------------------------------------------------------

# `Plots.Subplot` and `Plots.Axis` support `getindex` but not `haskey`, and throw for keys that are
# not in their defaults — hence the catch rather than a lookup.
function _get(d, k, default = nothing)
    v = try
        d[k]
    catch
        return default
    end
    v === nothing ? default : v
end

_string(::Nothing) = ""
_string(::Missing) = ""
_string(x) = string(x)
_string(x::Symbol) = String(x)

# Units and limits arrive as Unitful quantities; everything downstream wants bare Float64.
_bare(x::Unitful.Quantity) = Float64(ustrip(x))
_bare(x::Real) = Float64(x)
_bare(x) = nothing

_unit_of(a) = (u = _get(a, :unit); u === nothing ? "" : string(u))

function _lims(a)
    l = _get(a, :lims)
    (l isa Tuple && length(l) == 2) || return nothing
    lo, hi = _bare(l[1]), _bare(l[2])
    (lo === nothing || hi === nothing) ? nothing : Float64[lo, hi]
end

function _clims(sp)
    c = _get(sp, :clims)
    (c isa Tuple && length(c) == 2) || return nothing
    lo, hi = _bare(c[1]), _bare(c[2])
    (lo === nothing || hi === nothing) ? nothing : Float64[lo, hi]
end

# A colour attribute is either a single colour or a gradient. Both are reduced to a short, stable
# string: gradients to three sampled stops, which is enough to tell viridis from inferno from
# thermal without storing 256 RGBA values per series.
_hex2(x) = string(round(Int, clamp(x, 0, 1) * 255), base = 16, pad = 2)
_hex(c) = (rgba = Plots.plot_color(c); "#" * _hex2(rgba.r) * _hex2(rgba.g) * _hex2(rgba.b) *
                                       (rgba.alpha < 1 ? _hex2(rgba.alpha) : ""))

_color_string(::Nothing) = ""
_color_string(c::Symbol) = String(c)
_color_string(c::Integer) = "palette:$c"
function _color_string(c)
    try
        if c isa Plots.PlotUtils.ColorGradient
            return "grad(" * join([_hex(c[t]) for t in (0.0, 0.5, 1.0)], ",") * ")"
        end
        return _hex(c)
    catch
        return string(c)
    end
end

# Plots replaces an unset label with an auto-generated "y1", "y2", … . Those are not real legend
# entries, and Makie will never produce them, so they are recorded but excluded from
# `legend_entries`.
_is_auto_label(l::AbstractString) = occursin(r"^y\d+$", l)

# Series types that never appear in a legend, whatever their label.
const NON_LEGEND_TYPES = Set([:heatmap, :contour, :contour3d, :surface, :image])

# `Plots.Surface` wraps the matrix of a heatmap/contour/surface series.
_unwrap(v) = v isa Plots.Surface ? v.surf : v

function _n_connections(d)
    c = _get(d, :connections)
    c === nothing && return 0
    c isa AbstractVector && return length(c)
    c isa Tuple && return length(first(c))
    0
end

## ----------------------------------------------------------------------------------------------
## Plot object -> CaseSpec
## ----------------------------------------------------------------------------------------------

function _axis_spec(sp, legend_entries, colormap)
    AxisSpec(
        _string(_get(sp, :title, "")),
        _string(_get(sp[:xaxis], :guide, "")),
        _string(_get(sp[:yaxis], :guide, "")),
        _string(_get(sp[:zaxis], :guide, "")),
        _unit_of(sp[:xaxis]), _unit_of(sp[:yaxis]), _unit_of(sp[:zaxis]),
        _lims(sp[:xaxis]), _lims(sp[:yaxis]), _clims(sp),
        _string(_get(sp, :colorbar, "")),
        _string(_get(sp, :colorbar_title, "")),
        colormap,
        _string(_get(sp, :aspect_ratio, "")),
        _string(_get(sp, :projection, "")),
        _string(_get(sp, :legend_position, "")),
        legend_entries,
    )
end

function _series_spec(s, subplot_idx::Int)
    d = s.plotattributes
    SeriesSpec(
        subplot_idx,
        _string(_get(d, :seriestype, "")),
        _string(_get(d, :label, "")),
        _color_string(_get(d, :seriescolor)),
        _color_string(_get(d, :linecolor)),
        _bare(_get(d, :linewidth)),
        _bare(something(_get(d, :seriesalpha), _get(d, :linealpha), missing)),
        _bare(_get(d, :fillalpha)),
        _string(_get(d, :markershape, "")),
        _bare(_get(d, :markersize)),
        _n_connections(d),
        data_spec(_unwrap(_get(d, :x))),
        data_spec(_unwrap(_get(d, :y))),
        data_spec(_unwrap(_get(d, :z))),
    )
end

"""
    plots_case_spec(plt::Plots.Plot; id, group, call) -> CaseSpec

Introspect a Plots plot object into the golden spec of one case.
"""
function plots_case_spec(plt::Plots.Plot; id::AbstractString, group::AbstractString,
                         call::AbstractString, note::AbstractString = "")
    subplots = plt.subplots
    subplot_index = IdDict(sp => i for (i, sp) in enumerate(subplots))

    series = SeriesSpec[]
    for s in plt.series_list
        idx = get(subplot_index, s.plotattributes[:subplot], 1)
        push!(series, _series_spec(s, idx))
    end

    entries = [String[] for _ in subplots]
    for spec in series
        st = Symbol(spec.seriestype)
        (st in NON_LEGEND_TYPES || isempty(spec.label) || _is_auto_label(spec.label)) && continue
        push!(entries[spec.subplot], spec.label)
    end

    # Plots keeps the colour gradient on the series, not on the subplot (`sp[:seriescolor]` is
    # unset), so the axis-level colormap — the thing a `Colorbar` needs on the Makie side — is taken
    # from the first gradient-coloured series of that subplot.
    colormaps = ["" for _ in subplots]
    for spec in series
        isempty(colormaps[spec.subplot]) || continue
        startswith(spec.color, "grad(") && (colormaps[spec.subplot] = spec.color)
    end

    axes = AxisSpec[_axis_spec(sp, entries[i], colormaps[i]) for (i, sp) in enumerate(subplots)]
    CaseSpec(; id, group, call, backend = "plots", axes, series, note)
end

"""
    error_spec(id, group, call, err; note = "") -> CaseSpec

Spec recording that a case could not be plotted, so one broken case does not abort a whole run.
"""
error_spec(id, group, call, err; note = "") =
    CaseSpec(; id, group, call, backend = "plots", error = sprint(showerror, err), note)

end # module
