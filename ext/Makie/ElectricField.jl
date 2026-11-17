# Makie recipe for the electric field magnitude.
#
# This is the first recipe of the Plots.jl -> Makie.jl migration and serves as the template the
# remaining recipes are modelled on. It mirrors the behaviour of the Plots.jl recipe in
# `src/PlotRecipes/ElectricField.jl` (magnitude `norm.(ef)`) and the cross-section slicing that the
# Plots recipe delegated to `src/PlotRecipes/Potentials.jl`.
#
# The shared, backend-agnostic data-prep helpers (`get_crosssection_idx_and_value`,
# `get_extended_ticks`) live in the main package (`src/PlotRecipes/plot_utils.jl`,
# `src/Axes/DiscreteAxis.jl`) and are reused here.

# Field-line plotting is not ported yet; the Plots.jl recipe `plot_electric_fieldlines` remains the
# way to draw them. No Makie entry point is declared until it lands, so that the package never
# exports a name that has no methods.

"""
    electricfieldplot(ef::ElectricField; kwargs...)
    electricfieldplot!(ax, ef::ElectricField; kwargs...)

Plot the magnitude of the electric field `ef` on a cross section of the detector.

Exactly one cross-section selector must be given: `r`, `φ` or `z` for a cylindrical field,
`x`, `y` or `z` for a Cartesian one. The value is snapped to the closest grid tick.
"""
@recipe ElectricFieldPlot (ef,) begin
    "Colormap of the field magnitude."
    colormap = SSD_COLORMAPS.field
    "Limits of the color scale; `Makie.automatic` takes them from the data."
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
    "Overlay contour lines of the field magnitude."
    contours_equal_potential = false
    "Number of contour levels drawn by `contours_equal_potential` (Plots' `levels` default)."
    contour_levels = 15
    "Colormap of the contour overlay."
    contour_colormap = SSD_COLORMAPS.contours
    "Mirror a cylindrical `φ` cross section to show the full detector."
    full_det = false
    "Draw the colour-mapped field as an embedded bitmap at this scale when the figure is saved as
     SVG or PDF. A heatmap *is* a raster, and CairoMakie otherwise writes one vector path per cell:
     a refined ICPC potential comes to 400 MB of SVG. At 4x the result is indistinguishable from the
     vector one; `false` keeps it vector."
    rasterize = 4
end

# Make the generic entry points (`plot(ef)`, `plot!(ax, ef)`) route to this recipe — that is what
# keeps the existing `plot(sim.electric_field, φ = 30)` call sites in the docs working. The
# `convert_arguments` method is the identity: the recipe consumes the `ElectricField` itself, and
# without it Makie would try (and fail) to convert the field into point data.
Makie.plottype(::ElectricField{<:Any, 3}) = ElectricFieldPlot
Makie.convert_arguments(::Type{<:ElectricFieldPlot}, ef::ElectricField{<:Any, 3}) = (ef,)

# Pick a PolarAxis for the cylindrical z cross-section (φ–r plane), a regular Axis otherwise.
function Makie.preferred_axis_type(p::ElectricFieldPlot{<:Tuple{<:ElectricField{T, 3, Cylindrical}}}) where {T}
    cross_section, = get_crosssection_idx_and_value(p[1][].grid,
        _selector(p, :r), _selector(p, :φ), _selector(p, :z))
    cross_section == :z ? Makie.PolarAxis : Makie.Axis
end
Makie.preferred_axis_type(::ElectricFieldPlot{<:Tuple{<:ElectricField{T, 3, Cartesian}}}) where {T} = Makie.Axis

## ----------------------------------------------------------------------------------------------
## Shared data preparation
## ----------------------------------------------------------------------------------------------

# The magnitudes (V/m) the recipe draws, as a plain matrix.
#
# Both layers need it. The recipe draws it; the figure-level wrapper needs its extrema, because the
# `Colorbar` has to be given a colour range explicitly (plan §1.4) and the field magnitude is the one
# scalar field for which `_get_potential_plot_information` has no entry: the Plots recipe sets no
# `clims`, so GR scaled each panel to the slice it drew. Taking the whole field's extrema instead
# would label the bar with a range the picture never uses — on the Cartesian fixture the `x` cross
# section spans 0…5.1e5 V/m of a field whose maximum is 1.5e6, a factor of three.
#
# `ef.data` rather than `ef`: `ElectricField <: AbstractArray{T, 3}` declares the *scalar* type while
# its elements are static vectors, so a non-scalar `getindex` on it would try to build an `Array{T}`.
function _efield_slice(ef::ElectricField{T, 3, Cylindrical}, cross_section::Symbol, idx::Int,
                       value, full_det::Bool) where {T <: SSDFloat}
    if cross_section == :φ
        slice = norm.(ef.data[:, idx, :])
        if full_det
            _, idx_mirror, _ = get_crosssection_idx_and_value(ef.grid, missing, value + 180, missing)
            slice = cat(norm.(ef.data[end:-1:2, idx_mirror, :]), slice, dims = 1)
        end
        slice
    elseif cross_section == :r
        norm.(ef.data[idx, :, :])
    else # :z, in storage order (r, φ); the recipe transposes it for the PolarAxis
        norm.(ef.data[:, :, idx])
    end
end

function _efield_slice(ef::ElectricField{T, 3, Cartesian}, cross_section::Symbol, idx::Int,
                       _value, _full_det::Bool = false) where {T <: SSDFloat}
    cross_section == :x ? norm.(ef.data[idx, :, :]) :
    cross_section == :y ? norm.(ef.data[:, idx, :]) :
                          norm.(ef.data[:, :, idx])
end

# `Makie.automatic` on a heatmap resolves to these same limits, but only inside the backend and only
# at render time: `plot.colorrange[]` stays `Automatic()`, so a colorbar built next to it cannot read
# the scale off the plot, and the tier-1 gate records "no colorrange" instead of the one the picture
# uses. Resolving it here through `color_limits` — which reproduces Makie's resolution, degenerate
# case included — makes the figure self-consistent and the scale introspectable, and changes no pixel.
_efield_colorrange(p, slice) =
    p.colorrange[] === Makie.automatic ? color_limits(slice) : p.colorrange[]


## ----------------------------------------------------------------------------------------------
## The recipe
## ----------------------------------------------------------------------------------------------

function Makie.plot!(p::ElectricFieldPlot{<:Tuple{<:ElectricField{T, 3, Cylindrical}}}) where {T <: SSDFloat}
    ef = p[1][]
    grid = ef.grid
    cross_section, idx, value = get_crosssection_idx_and_value(grid, p.r[], p.φ[], p.z[])
    slice = _efield_slice(ef, cross_section, idx, value, p.full_det[])
    crange = _efield_colorrange(p, slice)

    if cross_section == :φ
        rs, rcenters = _edges(grid.r, internal_length_unit), _centers(grid.r, internal_length_unit)
        if p.full_det[]
            rs, rcenters = _mirrored_edges(grid.r, internal_length_unit),
                           _mirrored_centers(grid.r, internal_length_unit)
        end
        zs = _edges(grid.z, internal_length_unit)
        heatmap!(p, rs, zs, slice; colormap = p.colormap, colorrange = crange,
                 rasterize = p.rasterize[])
        if p.contours_equal_potential[]
            contour!(p, rcenters, _centers(grid.z, internal_length_unit), slice;
                     colormap = p.contour_colormap, levels = p.contour_levels)
        end
    elseif cross_section == :r
        heatmap!(p, display_values(_angular_edges(grid.φ), internal_angle_unit),
                 _edges(grid.z, internal_length_unit), slice;
                 rasterize = p.rasterize[],
            colormap = p.colormap, colorrange = crange)
    else # cross_section == :z  -> plotted into a PolarAxis (φ, r)
        polar_cells!(p, grid.φ, _edges(grid.r, internal_length_unit), slice;
                     colormap = p.colormap, colorrange = crange, rasterize = p.rasterize[])
    end
    return p
end

function Makie.plot!(p::ElectricFieldPlot{<:Tuple{<:ElectricField{T, 3, Cartesian}}}) where {T <: SSDFloat}
    ef = p[1][]
    grid = ef.grid
    cross_section, idx, _ = get_crosssection_idx_and_value(grid, p.x[], p.y[], p.z[])
    slice = _efield_slice(ef, cross_section, idx, nothing)
    crange = _efield_colorrange(p, slice)

    ax1, ax2 = if cross_section == :x
        grid.y, grid.z
    elseif cross_section == :y
        grid.x, grid.z
    else # :z
        grid.x, grid.y
    end
    heatmap!(p, _edges(ax1, internal_length_unit), _edges(ax2, internal_length_unit), slice;
             rasterize = p.rasterize[],
             colormap = p.colormap, colorrange = crange)
    if p.contours_equal_potential[]
        contour!(p, _centers(ax1, internal_length_unit), _centers(ax2, internal_length_unit), slice;
                 colormap = p.contour_colormap, levels = p.contour_levels)
    end
    return p
end


## ----------------------------------------------------------------------------------------------
## Figure-level wrapper
## ----------------------------------------------------------------------------------------------

# Same string as the Plots recipe, from the same helper, with the `sigdigits = 2` that recipe has
# always used — and in display units, so the title agrees with the axis labels (plan §2.2).
_efield_title(cross_section, value, units) =
    potential_plot_title(EFIELD_PLOT_NAME, cross_section, value, display_unit(units), sigdigits = 2)

function _efield_axis_setup(ef::ElectricField{T, 3, Cylindrical};
                            r = missing, φ = missing, z = missing,
                            full_det::Bool = false, kwargs...) where {T}
    grid = ef.grid
    cross_section, idx, value, units = get_crosssection_idx_and_value(grid, r, φ, z)
    common = (; title = _efield_title(cross_section, value, units),
              colorbar_label = unit_label(EFIELD_COLORBAR_NAME, internal_efield_unit; format = :parens),
              colorrange = color_limits(_efield_slice(ef, cross_section, idx, value, full_det)))
    lu = internal_length_unit
    if cross_section === :φ
        rlims = full_det ? (-last(grid.r), last(grid.r)) : extrema(grid.r)
        (; common..., axis = (; axistype = Makie.Axis, xname = "r", yname = "z",
           xunit = lu, yunit = lu, xlims = rlims, ylims = extrema(grid.z),
           aspect_ratio_one = true))
    elseif cross_section === :r
        (; common..., axis = (; axistype = Makie.Axis, xname = "φ", yname = "z",
           xunit = internal_angle_unit, yunit = lu, ylims = extrema(grid.z),
           xlims = extrema(_angular_edges(grid.φ)), xticks = _angle_ticks(grid.φ),
           _shell_aspect(value, extrema(grid.z))...))
    else
        under_labels = _efield_slice(ef, cross_section, idx, value, full_det)[:, 1]
        (; common..., axis = (; axistype = Makie.PolarAxis,
           rticklabelcolor = readable_on(under_labels, SSD_COLORMAPS.field, common.colorrange)))
    end
end

function _efield_axis_setup(ef::ElectricField{T, 3, Cartesian};
                            x = missing, y = missing, z = missing, kwargs...) where {T}
    grid = ef.grid
    cross_section, idx, value, units = get_crosssection_idx_and_value(grid, x, y, z)
    ax1, ax2, xname, yname = if cross_section === :x
        grid.y, grid.z, "y", "z"
    elseif cross_section === :y
        grid.x, grid.z, "x", "z"
    else
        grid.x, grid.y, "x", "y"
    end
    (; title = _efield_title(cross_section, value, units),
       colorbar_label = unit_label(EFIELD_COLORBAR_NAME, internal_efield_unit; format = :parens),
       colorrange = color_limits(_efield_slice(ef, cross_section, idx, value)),
       axis = (; axistype = Makie.Axis, xname, yname,
               xunit = internal_length_unit, yunit = internal_length_unit,
               xlims = extrema(ax1), ylims = extrema(ax2), aspect_ratio_one = true))
end

function electricfieldplot(gp::Makie.GridPosition, ef::ElectricField{<:Any, 3};
                           axis = (;), colorbar = (;), kwargs...)
    setup = _efield_axis_setup(ef; kwargs...)
    gl = panel_layout(gp)
    ax = ssd_axis(gl[1, 1]; title = setup.title, setup.axis..., axis...)
    plt = electricfieldplot!(ax, ef; kwargs...)
    crange = get(kwargs, :colorrange, Makie.automatic)
    crange === Makie.automatic && (crange = setup.colorrange)
    add_colorbar!(gl, plt; axis = ax, label = setup.colorbar_label, colorrange = crange,
                  colorbar...)
    Makie.AxisPlot(ax, plt)
end

function electricfieldplot(ef::ElectricField{<:Any, 3}; figure = (;), kwargs...)
    fig = Makie.Figure(; figure...)
    ap = electricfieldplot(fig[1, 1], ef; kwargs...)
    Makie.FigureAxisPlot(fig, ap.axis, ap.plot)
end

Makie.plot(ef::ElectricField{<:Any, 3}; kwargs...) = electricfieldplot(ef; kwargs...)


## ----------------------------------------------------------------------------------------------
## Electric field lines
##
## 12 of the ~142 parity cases. The computation - which cross section, where charges spawn, where
## they drift, how a path projects onto the plane - lives in `src/PlotRecipes/plot_utils.jl`, so both
## backends trace the same lines. `streamplot` is a phase-2 replacement for the algorithm, not a port
## of it.
##
## The Plots recipe emits one series per line, so a finely sampled detector produces hundreds. Here
## every line is one `lines!` with a non-finite vertex between them, the same choice the drift-path
## layer makes.
## ----------------------------------------------------------------------------------------------

"""
    fieldlinesplot(sim; kwargs...)
    fieldlinesplot!(ax, sim; kwargs...)

Draw the electric field lines of `sim` on one cross section.

Charges are spawned just inside the semiconductor at the contact surfaces and drifted parallel to
**E**, in both directions, so each line runs from contact to contact. `sim` must already carry an
electric field.

The lines are white, to be drawn over a field-magnitude plot:

```julia
fig = electricfieldplot(sim.electric_field; φ = 0, full_det = true)
fieldlinesplot!(fig.axis, sim; φ = 0, full_det = true)
```

The non-mutating form has no field underneath and draws them black. `color` overrides either.

This is the Makie counterpart of `plot_electric_fieldlines`, under a different name so that both
backends stay usable in one session.

## Keywords
* `r`, `φ`, `z` / `x`, `y`, `z`: the cross section. At most one; none selects `φ = 0`, or the middle
  `x` slice of a Cartesian grid.
* `sampling`: spacing of the surface sampling that places the charges.
* `offset`: how far inside the surface a charge starts. A charge left in the contact will not drift.
* `spacing`: keep every `n`-th line, for when fine sampling clutters the picture.
* `skip_contact`: the contact *not* sampled, so lines start on one polarity only.
* `max_nsteps`, `Δt`: the drift integration. Lines ending mid-detector want more steps or a smaller
  step size.
* `full_det`: annotate the title with the angle mirrored onto the left half.
* `color`, `linewidth`: the lines.
"""
@recipe FieldLinesPlot (simulation,) begin
    "Radial cross section (cylindrical grids)."
    r = missing
    "Angular cross section (cylindrical grids)."
    φ = missing
    "Cross section along x (Cartesian grids)."
    x = missing
    "Cross section along y (Cartesian grids)."
    y = missing
    "Cross section along z."
    z = missing
    "Maximum number of drift steps per line."
    max_nsteps = 5000
    "Spacing of the contact-surface sampling that places the charges."
    sampling = 2u"mm"
    "How far inside the surface a charge starts."
    offset = 0.5u"mm"
    "Keep every `n`-th line."
    spacing = 1
    "Drift step size."
    Δt = 1u"ns"
    "Annotate the title with the mirrored angle."
    full_det = false
    "The contact charges are *not* spawned on."
    skip_contact = 1
    "Colour of the lines."
    color = Makie.automatic
    "Width of the lines."
    linewidth = Makie.automatic
end

Makie.convert_arguments(::Type{<:FieldLinesPlot}, sim::Simulation) = (sim,)

_fieldline_selectors(p) = (; x = p.x[], y = p.y[], z = p.z[], r = p.r[], φ = p.φ[])
_fieldline_drift_options(p) = (; max_nsteps = p.max_nsteps[], sampling = p.sampling[],
                               offset = p.offset[], spacing = p.spacing[], Δt = p.Δt[],
                               skip_contact = p.skip_contact[])

# Every line end to end, separated by a non-finite vertex so Makie lifts the pen between them.
#
# The horizontal axis is an angle for the two cylindrical cross sections that have one, and there the
# two differ: a `PolarAxis` takes its θ in radians (the one exception in D7), while an `r =` cross
# section is a plain axis whose angle is shown in degrees like every other angle here.
function _fieldline_vertices(paths, dim_symbol::Symbol, value::Real, S)
    polar = S == Cylindrical && dim_symbol === :z
    xunit = dim_symbol === :r ? internal_angle_unit : internal_length_unit
    xs, ys = Float64[], Float64[]
    for path in paths
        isempty(path) && continue
        u, w = fieldline_projection(path, dim_symbol, value, S)
        if !isempty(xs)
            push!(xs, NaN)
            push!(ys, NaN)
        end
        append!(xs, polar ? u : display_values(u, xunit))
        append!(ys, display_values(w, internal_length_unit))
    end
    xs, ys
end

function Makie.plot!(p::FieldLinesPlot)
    sim = p[1][]
    S = get_coordinate_system(sim.electric_field.grid)
    dim_symbol, value, _ = fieldline_cross_section(sim; _fieldline_selectors(p)...)
    paths = fieldline_paths(sim, dim_symbol, value; _fieldline_drift_options(p)...)
    xs, ys = _fieldline_vertices(paths, dim_symbol, value, S)
    # A cross section that no contact surface meets spawns no charge and draws nothing. That is a
    # property of the geometry, not a failure: the CGD contacts are slabs perpendicular to x, so an
    # x cross section of it is legitimately empty.
    isempty(xs) && return p
    Makie.lines!(p, xs, ys;
                 color = color_or_default(p.color[], SSD_COLORS.fieldline),
                 linewidth = color_or_default(p.linewidth[], SSD_LINEWIDTHS.fieldline))
    return p
end


# Everything the wrapper has to set and the recipe cannot: which axis type, the labels, the units and
# the title. The cross sections match the field-magnitude ones, except that a `φ` cut is always drawn
# across the full diameter - the sample lines run from -r to +r, so a field line reaches both halves
# whatever `full_det` says, and `full_det` only annotates the title.
function _fieldlines_axis_setup(sim::Simulation; x = missing, y = missing, z = missing,
                               r = missing, φ = missing, full_det::Bool = false, kwargs...)
    grid = sim.electric_field.grid
    S = get_coordinate_system(grid)
    dim_symbol, value, units = fieldline_cross_section(sim; x, y, z, r, φ)
    title = fieldline_title(dim_symbol, value, display_unit(units);
                            full_det = fieldline_full_det(dim_symbol, full_det))
    lu = internal_length_unit
    axis = if S == Cylindrical
        if dim_symbol === :φ
            (; axistype = Makie.Axis, xname = "r", yname = "z", xunit = lu, yunit = lu,
               xlims = (-last(grid.r), last(grid.r)), ylims = extrema(grid.z),
               aspect_ratio_one = true)
        elseif dim_symbol === :r
            (; axistype = Makie.Axis, xname = "φ", yname = "z",
               xunit = internal_angle_unit, yunit = lu, ylims = extrema(grid.z),
               xlims = extrema(_angular_edges(grid.φ)), xticks = _angle_ticks(grid.φ),
           _shell_aspect(value, extrema(grid.z))...)
        else
            (; axistype = Makie.PolarAxis)
        end
    else
        # A cross section is drawn in the two dimensions it did *not* cut. The Plots recipe labels
        # the vertical axis "z" whatever the cut, which for a z cross section names the dimension it
        # just removed (plan F5); the labels here follow the data.
        ax1, ax2, xname, yname = dim_symbol === :x ? (grid.y, grid.z, "y", "z") :
                                 dim_symbol === :y ? (grid.x, grid.z, "x", "z") :
                                                     (grid.x, grid.y, "x", "y")
        (; axistype = Makie.Axis, xname, yname, xunit = lu, yunit = lu,
           xlims = extrema(ax1), ylims = extrema(ax2), aspect_ratio_one = true)
    end
    (; title, axis)
end

function fieldlinesplot(gp::Makie.GridPosition, sim::Simulation; axis = (;), kwargs...)
    setup = _fieldlines_axis_setup(sim; kwargs...)
    ax = ssd_axis(gp; title = setup.title, setup.axis..., axis...)
    # White is for lines over a field-magnitude plot. This axis is empty, so they are drawn dark.
    plt = fieldlinesplot!(ax, sim; color = SSD_COLORS.fieldline_alone, kwargs...)
    Makie.AxisPlot(ax, plt)
end

function fieldlinesplot(sim::Simulation; figure = (;), kwargs...)
    fig = Makie.Figure(; figure...)
    ap = fieldlinesplot(fig[1, 1], sim; kwargs...)
    Makie.FigureAxisPlot(fig, ap.axis, ap.plot)
end
