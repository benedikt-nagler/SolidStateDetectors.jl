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
    colormap = :inferno
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
    "Mirror a cylindrical `φ` cross section to show the full detector."
    full_det = false
end

# Make the generic entry points (`plot(ef)`, `plot!(ax, ef)`) route to this recipe — that is what
# keeps the existing `plot(sim.electric_field, φ = 30)` call sites in the docs working. The
# `convert_arguments` method is the identity: the recipe consumes the `ElectricField` itself, and
# without it Makie would try (and fail) to convert the field into point data.
Makie.plottype(::ElectricField{<:Any, 3}) = ElectricFieldPlot
Makie.convert_arguments(::Type{<:ElectricFieldPlot}, ef::ElectricField{<:Any, 3}) = (ef,)

# `preferred_axis_type` runs before the recipe's default attributes are merged into the plot, so the
# only place the user-passed cross-section keywords exist yet is the raw kwargs dict `p.kw`. An absent
# key means "not selected" (missing) — exactly the default used by `get_crosssection_idx_and_value`.
_selector(p, key) = Makie.to_value(get(p.kw, key, missing))

# Pick a PolarAxis for the cylindrical z cross-section (φ–r plane), a regular Axis otherwise.
function Makie.preferred_axis_type(p::ElectricFieldPlot{<:Tuple{<:ElectricField{T, 3, Cylindrical}}}) where {T}
    cross_section, = get_crosssection_idx_and_value(p[1][].grid,
        _selector(p, :r), _selector(p, :φ), _selector(p, :z))
    cross_section == :z ? Makie.PolarAxis : Makie.Axis
end
Makie.preferred_axis_type(::ElectricFieldPlot{<:Tuple{<:ElectricField{T, 3, Cartesian}}}) where {T} = Makie.Axis

function Makie.plot!(p::ElectricFieldPlot{<:Tuple{<:ElectricField{T, 3, Cylindrical}}}) where {T <: SSDFloat}
    ef = p[1][]
    grid = ef.grid
    ef_magn = norm.(ef) # plain magnitudes (V/m)
    cross_section, idx, value = get_crosssection_idx_and_value(grid, p.r[], p.φ[], p.z[])

    if cross_section == :φ
        if p.full_det[]
            _, idx_mirror, _ = get_crosssection_idx_and_value(grid, missing, value + 180, missing)
            data = cat(ef_magn[end:-1:2, idx_mirror, :], ef_magn[:, idx, :], dims = 1)
            xs = vcat(-1 .* reverse(grid.r[2:end]), collect(grid.r))
            heatmap!(p, xs, collect(grid.z), data; colormap = p.colormap, colorrange = p.colorrange)
        else
            heatmap!(p, midpoints(get_extended_ticks(grid.r)), midpoints(get_extended_ticks(grid.z)),
                ef_magn[:, idx, :]; colormap = p.colormap, colorrange = p.colorrange)
        end
        if p.contours_equal_potential[]
            contour!(p, collect(grid.r), collect(grid.z), ef_magn[:, idx, :]; colormap = :thermal)
        end
    elseif cross_section == :r
        heatmap!(p, collect(grid.φ), collect(grid.z), ef_magn[idx, :, :];
            colormap = p.colormap, colorrange = p.colorrange)
    else # cross_section == :z  -> plotted into a PolarAxis (φ, r)
        heatmap!(p, collect(grid.φ), collect(grid.r), ef_magn[:, :, idx]';
            colormap = p.colormap, colorrange = p.colorrange)
    end
    return p
end

function Makie.plot!(p::ElectricFieldPlot{<:Tuple{<:ElectricField{T, 3, Cartesian}}}) where {T <: SSDFloat}
    ef = p[1][]
    grid = ef.grid
    ef_magn = norm.(ef)
    cross_section, idx, _ = get_crosssection_idx_and_value(grid, p.x[], p.y[], p.z[])

    xs, ys, data = if cross_section == :x
        collect(grid.y), collect(grid.z), ef_magn[idx, :, :]
    elseif cross_section == :y
        collect(grid.x), collect(grid.z), ef_magn[:, idx, :]
    else # :z
        collect(grid.x), collect(grid.y), ef_magn[:, :, idx]
    end
    heatmap!(p, xs, ys, data; colormap = p.colormap, colorrange = p.colorrange)
    return p
end
