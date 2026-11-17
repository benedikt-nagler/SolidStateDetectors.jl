# Makie recipe and figure-level wrapper for the scalar fields on a simulation grid:
# `ElectricPotential`, `WeightingPotential`, `EffectiveChargeDensity`, `PointTypes`, `ImpurityScale`
# and `DielectricDistribution` — 49 of the ~142 parity cases, and the mirror image of the two Plots
# recipes in `src/PlotRecipes/Potentials.jl`.
#
# Two layers (plan §1.1 / D1), because a Makie recipe owns nothing outside itself:
#
#   * `ScalarPotentialPlot` (`scalarpotentialplot!`) draws the heatmap and the optional contour
#     overlay into whatever axis it is handed;
#   * `scalarpotentialplot(sp; …)` builds the `Figure`, the `Axis`/`PolarAxis` (title, labels,
#     limits, `DataAspect`) and the `Colorbar` — everything the Plots recipe set at subplot level.
#     `Makie.plot(sp; …)` routes here too, so the existing `plot(sim.electric_potential, φ = 0)`
#     call sites keep producing a titled plot with a colorbar rather than a bare heatmap.
#
# Three things that are deliberately *not* a literal translation of the Plots recipe:
#
#   * **No transposes.** Plots' `heatmap(x, y, z)` wants `size(z) == (length(y), length(x))`, Makie's
#     wants `(length(x), length(y))`, and SSD stores `[r, φ, z]` / `[x, y, z]` — so the `'` in every
#     Plots slice disappears here. The exception is the cylindrical `z` cross section, which is
#     `(r, φ)` in storage and `(φ, r)` on a `PolarAxis` (plan §1.7).
#   * **Cell edges, not centers.** Plots was handed `midpoints(midpoints(get_extended_ticks(ax)))`,
#     i.e. `n` centers, and clipped the outermost half-cells with `xlims`. Makie accepts `n+1`
#     *edges*, which is what `get_extended_ticks` produces one `midpoints` earlier, so cell
#     boundaries land exactly on the grid midpoints. The limits still clip the outer cells, which on
#     an `:infinite`/`:r0` axis extend far past the detector by construction.
#   * **Display units in the data.** Makie 0.24's `heatmap` cannot take Unitful data (see
#     `display_values` in `units.jl`), so lengths are rescaled to mm and angles to degrees here and
#     the unit goes into the axis label (D2). The one exception is the `θ` argument of a
#     `PolarAxis`, which Makie requires in radians.

"""
    scalarpotentialplot(sp::ScalarPotential; kwargs...)
    scalarpotentialplot!(ax, sp::ScalarPotential; kwargs...)

Plot a cross section of a scalar field of a [`Simulation`](@ref) — an `ElectricPotential`,
`WeightingPotential`, `EffectiveChargeDensity`, `PointTypes`, `ImpurityScale` or
`DielectricDistribution`.

Exactly one cross-section selector may be given: `r`, `φ` or `z` for a cylindrical grid, `x`, `y` or
`z` for a Cartesian one; the value is snapped to the closest grid tick, and the default is `φ = 0`
(cylindrical) / `x = 0` (Cartesian). A cylindrical `z` cross section is drawn into a `PolarAxis`.

The non-mutating form returns a `FigureAxisPlot` whose figure carries the title, the axis labels and
the colorbar; the mutating form draws into an existing axis and sets none of them.

## Keywords
* `contours_equal_potential`: overlay contour lines of the plotted field.
* `full_det`: mirror a cylindrical `φ` cross section about `r = 0` to show the full detector.
* `all_point_types`: for `PointTypes`, plot the full bit mask instead of only its lowest three bits.
* `colormap`, `colorrange`: default to the field's own gradient and limits (`:viridis` on
  `(0, 1)` for a weighting potential, and so on).
* `figure`, `axis`, `colorbar`: keywords forwarded to the `Figure`, the axis and the `Colorbar`.
"""
@recipe ScalarPotentialPlot (sp,) begin
    "Colormap; `Makie.automatic` takes the gradient the field is conventionally drawn with."
    colormap = Makie.automatic
    "Limits of the color scale; `Makie.automatic` takes the field's own limits."
    colorrange = Makie.automatic
    "Cylindrical cross section at this radius."
    r = missing
    "Cylindrical cross section at this angle (in degrees, unless a unit is attached)."
    φ = missing
    "Cross section at this height; cylindrical fields are drawn into a `PolarAxis`."
    z = missing
    "Cartesian cross section at this `x`."
    x = missing
    "Cartesian cross section at this `y`."
    y = missing
    "Overlay contour lines of equal potential."
    contours_equal_potential = false
    "Number of contour levels drawn by `contours_equal_potential` (Plots' `levels` default)."
    contour_levels = 15
    "Colormap of the contour overlay."
    contour_colormap = SSD_COLORMAPS.contours
    "Mirror a cylindrical `φ` cross section to show the full detector."
    full_det = false
    "Plot the full `PointTypes` bit mask rather than only the lowest three bits."
    all_point_types = false
    "Draw the colour-mapped field as an embedded bitmap at this scale when the figure is saved as
     SVG or PDF. A heatmap *is* a raster, and CairoMakie otherwise writes one vector path per cell:
     a refined ICPC potential comes to 400 MB of SVG. At 4x the result is indistinguishable from the
     vector one; `false` keeps it vector."
    rasterize = 4
end

# `plot(sp)` / `plot!(ax, sp)` must route to the recipe rather than to Makie's `AbstractArray`
# fallbacks — a `ScalarPotential` *is* an `AbstractArray{T, 3}`, so without these it would be
# interpreted as raw volume data.
Makie.plottype(::ScalarPotential{<:Any, 3}) = ScalarPotentialPlot
Makie.convert_arguments(::Type{<:ScalarPotentialPlot}, sp::ScalarPotential{<:Any, 3}) = (sp,)

function Makie.preferred_axis_type(p::ScalarPotentialPlot{<:Tuple{<:ScalarPotential{T, 3, Cylindrical}}}) where {T}
    cross_section, = get_crosssection_idx_and_value(p[1][].grid,
        _selector(p, :r), _selector(p, :φ), _selector(p, :z))
    cross_section == :z ? Makie.PolarAxis : Makie.Axis
end
Makie.preferred_axis_type(::ScalarPotentialPlot{<:Tuple{<:ScalarPotential{T, 3, Cartesian}}}) where {T} = Makie.Axis


## ----------------------------------------------------------------------------------------------
## Shared data preparation
## ----------------------------------------------------------------------------------------------

# A cylindrical field covering less than 2π is expanded, exactly as the Plots recipe does, so that a
# `z` cross section is a full disc and `full_det` has a mirror side to draw. `get_2π_potential` has
# no method for a `DielectricDistribution` (there is no `ScalarPotential(::DielectricDistribution,
# data, grid)` constructor), which is why the Plots recipe leaves that one alone too.
_full_φ(sp::ScalarPotential{T, 3, Cylindrical}) where {T} =
    (sp.grid[2][end] - sp.grid[2][1] ≈ 2π) ? sp : get_2π_potential(sp, n_points_in_φ = 72)
_full_φ(ϵ::DielectricDistribution{T, 3, Cylindrical}) where {T} = ϵ
_full_φ(sp::ScalarPotential{T, 3, Cartesian}) where {T} = sp

# The dielectric distribution lives on the *extended* grid: its first and last ticks are the solver's
# boundary points, which sit far outside the detector. The cylindrical Plots recipe therefore trims
# the limits to the inner ticks and drops the outer radial cells of a `z` cross section; the
# Cartesian one does neither, so its limits stay wide — kept as-is here for parity, and flagged for
# the tier-2 review as a candidate improvement.
_trim_radial(::ScalarPotential) = false
_trim_radial(::DielectricDistribution) = true

# Point types are a bit mask: by default only the lowest three bits (undepleted / depleted / fixed)
# are plotted, which is what makes the `(0, 7)` color range meaningful.
_field_values(sp::ScalarPotential, ::Bool) = sp.data
_field_values(pt::PointTypes, all_point_types::Bool) = all_point_types ? pt.data : pt.data .& 0x07

# `PointType` is a `UInt8`; Makie wants floats.
_float_values(a::AbstractArray{<:AbstractFloat}) = a
_float_values(a::AbstractArray) = Float32.(a)

# `color_limits` rather than a bare conversion: `_get_potential_plot_information` returns the field's
# own extrema for an `ElectricPotential`, an `EffectiveChargeDensity` or a `DielectricDistribution`,
# and those collapse to a single value on a constant field — where an explicit `colorrange` throws.
_clim_numbers(clims) = color_limits(Float64(ustrip(clims[1])), Float64(ustrip(clims[2])))

function _color_defaults(p, sp)
    gradient, clims, _, _ = _get_potential_plot_information(sp, p.all_point_types[])
    cmap = p.colormap[] === Makie.automatic ? gradient : p.colormap[]
    crange = p.colorrange[] === Makie.automatic ? _clim_numbers(clims) : p.colorrange[]
    cmap, crange
end

# `n+1` cell edges of a grid axis, in display units. `get_extended_ticks` adds the boundary ticks
# the SOR solver uses, and their midpoints are the cell boundaries of the `n` data points.
#
# The two outer edges are clamped to the axis' own tick range: on an `:infinite` or `:r0` axis
# `get_extended_ticks` extends by the *full* axis width, so the outermost cell would otherwise reach
# far outside the detector — half again as far as its radius. On a Cartesian axis the limits used to
# hide that; a `PolarAxis` has no limits to hide it with, and drew a detector 50% too wide.
function _edges(ax, u = nothing)
    ticks = collect(ax)
    display_values(clamp.(midpoints(get_extended_ticks(ax)), first(ticks), last(ticks)), u)
end

# Cell centers, in display units — what `contour!` needs, since a contour is drawn *through* the
# data points rather than around them.
_centers(ax, u = nothing) = display_values(collect(ax), u)

# The same two, mirrored about zero for `full_det`: `2n-1` centers and the `2n` edges around them.
# The outermost edges are the plain ones reflected, so the mirrored picture is exactly symmetric.
_mirrored_centers(ax, u = nothing) = (t = collect(ax); display_values(vcat(-reverse(t[2:end]), t), u))
function _mirrored_edges(ax, u = nothing)
    outer = last(collect(ax))
    display_values(vcat(-outer, midpoints(_mirrored_centers(ax)), outer), u)
end

# The angular cell boundaries of a grid axis, in radians: midway between neighbouring ticks. A grid
# that covers the full circle is periodic, so its first and last cell are as wide as the others and
# the boundaries run from half a cell below the first tick to half a cell above the last. A partial
# range ends at its own first and last tick, like every other axis.
function _angular_edges(φaxis)
    φ = Float64.(collect(φaxis))
    n = length(φ)
    n == 1 && return [φ[1] - π, φ[1] + π]
    step = (φ[end] - φ[1]) / (n - 1)
    periodic = isapprox(φ[end] - φ[1] + step, 2π; rtol = 1e-3)
    periodic || return clamp.(midpoints(get_extended_ticks(φaxis)), φ[1], φ[end])
    vcat(φ[1] - step / 2, midpoints(φ), φ[end] + step / 2)
end

# A polar cross section is drawn as one flat quad per grid cell, in (θ, r).
#
# Not `heatmap!`: its cells are axis-aligned rectangles, which the polar transform collapses into a
# sliver (Makie 0.24.13). Not `surface!` either: it interpolates colours between grid points, which
# smears a refined grid into a glow and invents point types that are not in the data.
#
# The edge of a quad is a chord, so an angular cell wider than `POLAR_MAX_CELL` is split into equal
# parts of the same value. A φ-symmetric field has one or two φ ticks and would otherwise be drawn as
# a line or a square.
const POLAR_MAX_CELL = 2π / 72

function polar_cells!(p, φaxis, r_edges::AbstractVector, slice::AbstractMatrix; kwargs...)
    θ_edges = _angular_edges(φaxis)
    pts = Makie.Point2f[]
    faces = Makie.GeometryBasics.GLTriangleFace[]
    values = Float32[]
    for j in axes(slice, 2)
        θ0, θ1 = θ_edges[j], θ_edges[j + 1]
        # The margin absorbs the rounding of `Float32` ticks, which would split a cell of exactly
        # the largest width in two.
        parts = max(1, ceil(Int, (θ1 - θ0) / POLAR_MAX_CELL - 0.01))
        for k in 1:parts
            a = θ0 + (θ1 - θ0) * (k - 1) / parts
            b = θ0 + (θ1 - θ0) * k / parts
            for i in axes(slice, 1)
                n = length(pts)
                push!(pts, Makie.Point2f(a, r_edges[i]), Makie.Point2f(b, r_edges[i]),
                      Makie.Point2f(b, r_edges[i + 1]), Makie.Point2f(a, r_edges[i + 1]))
                push!(faces, Makie.GeometryBasics.GLTriangleFace(n + 1, n + 2, n + 3),
                      Makie.GeometryBasics.GLTriangleFace(n + 1, n + 3, n + 4))
                v = Float32(slice[i, j])
                push!(values, v, v, v, v)
            end
        end
    end
    Makie.mesh!(p, pts, faces; color = values, shading = Makie.NoShading, kwargs...)
end


## ----------------------------------------------------------------------------------------------
## The recipe
## ----------------------------------------------------------------------------------------------

function Makie.plot!(p::ScalarPotentialPlot{<:Tuple{<:ScalarPotential{T, 3, Cylindrical}}}) where {T <: SSDFloat}
    sp = _full_φ(p[1][])
    grid = sp.grid
    cross_section, idx, value = get_crosssection_idx_and_value(grid, p.r[], p.φ[], p.z[])
    data = _float_values(_field_values(sp, p.all_point_types[]))
    cmap, crange = _color_defaults(p, sp)

    if cross_section == :φ
        slice = data[:, idx, :]
        rs, rcenters = _edges(grid.r, internal_length_unit), _centers(grid.r, internal_length_unit)
        if p.full_det[]
            # The mirror side is the cross section 180° away, with its radial axis reversed. Both
            # the heatmap and the contour overlay are mirrored — in Plots only because both `@series`
            # blocks repeat the same `cat`; here they share one array.
            _, idx_mirror, _ = get_crosssection_idx_and_value(grid, missing, value + 180, missing)
            slice = cat(data[end:-1:2, idx_mirror, :], slice, dims = 1)
            rs, rcenters = _mirrored_edges(grid.r, internal_length_unit),
                           _mirrored_centers(grid.r, internal_length_unit)
        end
        zs = _edges(grid.z, internal_length_unit)
        Makie.heatmap!(p, rs, zs, slice; colormap = cmap, colorrange = crange,
                       rasterize = p.rasterize[])
        if p.contours_equal_potential[]
            Makie.contour!(p, rcenters, _centers(grid.z, internal_length_unit), slice;
                           colormap = p.contour_colormap, levels = p.contour_levels)
        end
    elseif cross_section == :r
        Makie.heatmap!(p, display_values(_angular_edges(grid.φ), internal_angle_unit),
                       _edges(grid.z, internal_length_unit),
                       data[idx, :, :]; colormap = cmap, colorrange = crange,
                       rasterize = p.rasterize[])
    else # :z — a PolarAxis, whose θ argument stays in radians
        # The dielectric's extended grid puts its outermost radial cells outside the detector, and
        # the Plots recipe drops them (`ϵ.data[2:end-1, :, idx]`) rather than clipping — a PolarAxis
        # has no limits to clip with.
        rows = _trim_radial(sp) ? (2:(size(data, 1) - 1)) : (1:size(data, 1))
        r_edges = _edges(grid.r, internal_length_unit)[first(rows):(last(rows) + 1)]
        polar_cells!(p, grid.φ, r_edges, data[rows, :, idx];
                     colormap = cmap, colorrange = crange, rasterize = p.rasterize[])
    end
    return p
end

function Makie.plot!(p::ScalarPotentialPlot{<:Tuple{<:ScalarPotential{T, 3, Cartesian}}}) where {T <: SSDFloat}
    sp = p[1][]
    grid = sp.grid
    cross_section, idx, _ = get_crosssection_idx_and_value(grid, p.x[], p.y[], p.z[])
    data = _float_values(_field_values(sp, p.all_point_types[]))
    cmap, crange = _color_defaults(p, sp)

    ax1, ax2, slice = if cross_section == :x
        grid.y, grid.z, data[idx, :, :]
    elseif cross_section == :y
        grid.x, grid.z, data[:, idx, :]
    else # :z
        grid.x, grid.y, data[:, :, idx]
    end
    u = internal_length_unit
    Makie.heatmap!(p, _edges(ax1, u), _edges(ax2, u), slice; colormap = cmap,
                   colorrange = crange, rasterize = p.rasterize[])
    if p.contours_equal_potential[]
        Makie.contour!(p, _centers(ax1, u), _centers(ax2, u), slice;
                       colormap = p.contour_colormap, levels = p.contour_levels)
    end
    return p
end


## ----------------------------------------------------------------------------------------------
## Figure-level wrapper
## ----------------------------------------------------------------------------------------------

# Axis limits of a cylindrical cross section: trimmed to the inner ticks for the dielectric
# distribution, whose grid carries the solver's boundary points (see `_trim_radial`).
_cyl_limits(ax, sp) = _trim_radial(sp) ? (ax[2], ax[end - 1]) : extrema(ax)

# Everything the wrapper has to know and the recipe cannot express: which axis type, what the axis is
# called, in which unit, over which range, and how the color scale is labelled.
function _axis_setup(sp::ScalarPotential{T, 3, Cylindrical}; r = missing, φ = missing, z = missing,
                     full_det::Bool = false, all_point_types::Bool = false, kwargs...) where {T}
    grid = sp.grid
    cross_section, _, value, units = get_crosssection_idx_and_value(grid, r, φ, z)
    _, clims, name, punit = _get_potential_plot_information(sp, all_point_types)
    common = (; title = potential_plot_title(name, cross_section, value, display_unit(units)),
              colorbar_label = unit_label(name, punit; format = :parens),
              colorrange = _clim_numbers(clims))
    lu = internal_length_unit
    if cross_section == :φ
        rlims = full_det ? (-last(grid.r), last(grid.r)) : _cyl_limits(grid.r, sp)
        (; common..., axis = (; axistype = Makie.Axis, xname = "r", yname = "z",
           xunit = lu, yunit = lu, xlims = rlims, ylims = _cyl_limits(grid.z, sp),
           aspect_ratio_one = true))
    elseif cross_section == :r
        (; common..., axis = (; axistype = Makie.Axis, xname = "φ", yname = "z",
           xunit = internal_angle_unit, yunit = lu, ylims = _cyl_limits(grid.z, sp),
           xlims = extrema(_angular_edges(grid.φ)), xticks = _angle_ticks(grid.φ),
           _shell_aspect(value, _cyl_limits(grid.z, sp))...))
    else # :z — a PolarAxis has no labels, no aspect and no Cartesian limits
        gradient, _, _, _ = _get_potential_plot_information(sp, all_point_types)
        _, idx, _ = get_crosssection_idx_and_value(grid, r, φ, z)
        under_labels = _float_values(_field_values(sp, all_point_types))[:, 1, idx]
        (; common..., axis = (; axistype = Makie.PolarAxis,
           rticklabelcolor = readable_on(under_labels, gradient, common.colorrange)))
    end
end

# An `r =` cross section is not a cut through the detector: it is the cylindrical *shell* at that
# radius, unrolled into (φ, z). Degrees and millimetres are not commensurable, so `DataAspect()`
# would be meaningless — but the shell has real proportions, `2πr` wide by `Δz` tall, and giving the
# axis box those proportions is what stops a 40 mm detector from being drawn as a full-width band.
# The ticks stay in degrees, which is what locates a contact.
# Quadrant ticks for the angular axis. A thin shell gets a narrow box, and Makie's automatic ticks
# then overlap; multiples of 90° are both readable at any width and what one actually looks for on a
# segmented detector. Derived from the grid's own φ range, which need not start at 0.
function _angle_ticks(φaxis, step = 90)
    lo, hi = extrema(display_values(collect(φaxis), internal_angle_unit))
    ticks = collect((step * ceil(lo / step)):step:(step * floor(hi / step)))
    length(ticks) >= 2 ? ticks : Makie.automatic
end

function _shell_aspect(r, zlims)
    circumference = 2π * display_values(r, internal_length_unit)
    height = display_values(zlims[2], internal_length_unit) -
             display_values(zlims[1], internal_length_unit)
    (circumference > 0 && height > 0) || return NamedTuple()
    (; aspect = Makie.AxisAspect(circumference / height))
end

function _axis_setup(sp::ScalarPotential{T, 3, Cartesian}; x = missing, y = missing, z = missing,
                     all_point_types::Bool = false, kwargs...) where {T}
    grid = sp.grid
    cross_section, _, value, units = get_crosssection_idx_and_value(grid, x, y, z)
    _, clims, name, punit = _get_potential_plot_information(sp, all_point_types)
    ax1, ax2, xname, yname = if cross_section == :x
        grid.y, grid.z, "y", "z"
    elseif cross_section == :y
        grid.x, grid.z, "x", "z"
    else
        grid.x, grid.y, "x", "y"
    end
    lu = internal_length_unit
    (; title = potential_plot_title(name, cross_section, value, display_unit(units)),
       colorbar_label = unit_label(name, punit; format = :parens),
       colorrange = _clim_numbers(clims),
       axis = (; axistype = Makie.Axis, xname, yname, xunit = lu, yunit = lu,
               xlims = extrema(ax1), ylims = extrema(ax2),
               aspect_ratio_one = true))
end

# A method on the recipe's own function: `(sp::ScalarPotential{<:Any, 3})` is more specific than the
# `(args...)` method `@recipe` generates, so this is the one a caller of `scalarpotentialplot(sp)`
# gets — no shadowing, no second function object.
function scalarpotentialplot(gp::Makie.GridPosition, sp::ScalarPotential{<:Any, 3};
                             axis = (;), colorbar = (;), kwargs...)
    sp_full = _full_φ(sp)
    setup = _axis_setup(sp_full; kwargs...)
    gl = panel_layout(gp)
    ax = ssd_axis(gl[1, 1]; title = setup.title, setup.axis..., axis...)
    plt = scalarpotentialplot!(ax, sp_full; kwargs...)
    # The colorbar has to be told the range explicitly (plan §1.4), so a caller-supplied
    # `colorrange` has to be picked up here as well or the bar would describe a different scale
    # from the plot it sits next to.
    crange = get(kwargs, :colorrange, Makie.automatic)
    crange === Makie.automatic && (crange = setup.colorrange)
    add_colorbar!(gl, plt; axis = ax, label = setup.colorbar_label, colorrange = crange,
                  colorbar...)
    Makie.AxisPlot(ax, plt)
end

function scalarpotentialplot(sp::ScalarPotential{<:Any, 3}; figure = (;), kwargs...)
    fig = Makie.Figure(; figure...)
    ap = scalarpotentialplot(fig[1, 1], sp; kwargs...)
    Makie.FigureAxisPlot(fig, ap.axis, ap.plot)
end

# `plot(sp)` is what the docs and existing user code call, and it should keep producing the titled,
# colorbar-carrying figure Plots produced — so it goes to the wrapper, not to the bare recipe.
# `plot!(ax, sp)` stays on the recipe, since an existing axis belongs to the caller.
Makie.plot(sp::ScalarPotential{<:Any, 3}; kwargs...) = scalarpotentialplot(sp; kwargs...)
