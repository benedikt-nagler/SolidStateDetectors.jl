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

# Field-line plotting is not yet ported to Makie (the `ElectricFieldLines` stub is still pending).
# Until it lands, the Plots.jl recipe `plot_electric_fieldlines` remains available. Defining these
# here gives a clear message instead of an `UndefVarError` from the main-module forwarding stub.
const _fieldlines_todo = "`electricfieldlinesplot` is not yet implemented in the Makie extension; " *
    "use the Plots.jl recipe `plot_electric_fieldlines` for now."
electricfieldlinesplot(args...; kwargs...) = error(_fieldlines_todo)
electricfieldlinesplot!(args...; kwargs...) = error(_fieldlines_todo)

@recipe(ElectricFieldPlot) do scene
    Attributes(
        colormap = :inferno,
        colorrange = Makie.automatic,
        # cross-section selectors (only one may be set); cylindrical: r/φ/z, Cartesian: x/y/z
        r = missing, φ = missing, z = missing, x = missing, y = missing,
        contours_equal_potential = false,
        full_det = false,
    )
end

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
