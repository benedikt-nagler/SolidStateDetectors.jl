# Backend-agnostic plotting helpers.
#
# These functions contain no RecipesBase/Plots (or Makie) code: they only do data preparation
# (cross-section lookup, plot metadata, field-line sampling) that is shared by both the existing
# Plots.jl recipes and the Makie extension. They live here (rather than inside the `@recipe` files)
# so they survive the eventual removal of the Plots.jl recipes.

## Legend labels
#
# Plots builds its legend from every series' `label`; Makie collects entries only from direct
# children of the axis, so a recipe that draws all contacts internally produces none and the
# figure-level wrapper has to build the `Legend` from explicit elements. Both backends must
# nevertheless print the *same* text, so the strings live here rather than in either recipe.

"""
    plot_label(obj) -> String

The legend label of a detector part: `"Semiconductor"`, `"<name> (id: <id>)"` or
`"Contact - id: <id>"` for a [`Contact`](@ref), the passive's name or `"Passive <id>"` for a
[`Passive`](@ref).
"""
plot_label(::Semiconductor) = "Semiconductor"
plot_label(c::Contact) = c.name != "" ? "$(c.name) (id: $(c.id))" : "Contact - id: $(c.id)"
plot_label(p::Passive) = p.name != "" ? p.name : "Passive $(p.id)"

# A charge cloud is labelled by the shell shape it is built from, a whole `NBodyChargeCloud` by its
# own name (its shells carry no label of their own), and a particle source by its type.
plot_label(::PlatonicSolid{N}) where {N} = PLATONIC_SOLID_NAMES[N]
plot_label(::NBodyChargeCloud) = "NBodyChargeCloud"
plot_label(m::AbstractParticleSource) = string(typeof(m).name.name)

# The platonic solids are all aliases of `PlatonicSolid{N}`, so the vertex count is the only thing
# that distinguishes a `Tetrahedron` from a `Hexahedron` — including for the name.
const PLATONIC_SOLID_NAMES = Dict{Int, String}(
    1 => "PointCharge", 4 => "Tetrahedron", 6 => "Octahedron",
    8 => "Hexahedron", 12 => "Icosahedron", 20 => "Dodecahedron",
)

# Labels of the two charge-carrier species of an `EHDriftPath`.
const electron_label = "e-"
const hole_label = "h+"

# Labels of the three curves of an electron/hole contribution plot, and the panel title of a grid
# plot. Both backends print these, so they live here for the same reason `plot_label` does.
const electron_contribution_label = "Electron contribution"
const hole_contribution_label = "Hole contribution"
const sum_contribution_label = "Sum"

"""
    grid_panel_title(dim) -> String

The heading of the grid-point-density panel for axis `dim` of a [`Grid`](@ref).

The Plots recipe writes it into the *x* guide, where it sits over ticks that are really a coordinate;
the Makie wrapper makes it the axis title and labels the coordinate instead.
"""
grid_panel_title(dim::Integer) = "Grid point density - Axis $(dim)"

"""
    waveform_signal_label(signal_unit) -> String

The name of the quantity a waveform's samples measure: `"Charge"` when they are calibrated in a
charge unit, `"Energy"` when in an energy unit, and `"signal"` otherwise.
"""
function waveform_signal_label(u::Unitful.Units)
    Unitful.dimension(u) === Unitful.dimension(internal_charge_unit) && return "Charge"
    Unitful.dimension(u) === Unitful.dimension(internal_energy_unit) && return "Energy"
    "signal"
end


## Charge-cloud connectivity

"""
    charge_cloud_connections(solid) -> Vector{Vector{Int}}

The edges of a [`PlatonicSolid`](@ref), as index paths into its `locations`.

Each entry is one polyline through the vertices, chosen so that the whole wireframe is drawn in as
few strokes as possible: `[1, 2, 3, 4, 2]` is a four-segment stroke, not four separate edges. A
`PointCharge` has no edges.

Shared by both plotting backends, so the two draw the same solid.

The method taking the type is for an [`NBodyChargeCloud`](@ref), whose `shell_structure` field holds
`PlatonicSolid{N, T}` itself rather than an instance of it.
"""
charge_cloud_connections(::PlatonicSolid{N}) where {N} = PLATONIC_SOLID_CONNECTIONS[N]
charge_cloud_connections(::Type{<:PlatonicSolid{N}}) where {N} = PLATONIC_SOLID_CONNECTIONS[N]

const PLATONIC_SOLID_CONNECTIONS = Dict{Int, Vector{Vector{Int}}}(
    1  => Vector{Int}[],
    4  => [[1, 2, 3, 4, 2], [3, 1, 4]],
    6  => [[1, 2, 3, 4, 5, 2, 6, 5, 1, 3, 6, 4, 1]],
    8  => [[1, 2, 5, 3, 6, 4, 7, 8, 5], [3, 1, 4], [2, 7], [6, 8]],
    12 => [[1, 2, 7, 8, 9, 10, 11, 7, 3, 8, 4, 9, 5, 10, 6, 11, 2, 3, 4, 5, 6, 2],
           [3, 1, 6], [5, 1, 4], [7, 12, 8], [9, 12, 10], [11, 12]],
    20 => [[1, 2, 3, 4, 5, 1, 6, 11, 7, 12, 8, 13, 9, 14, 10, 15, 6],
           [15, 20, 16, 17, 18, 19, 20],
           ([j, j + 5] for j in (2, 3, 4, 5, 11, 12, 13, 14))...],
)

"""
    charge_cloud_shell_markersize(markersize, shell) -> Float64

The marker size of the `shell`-th shell of an [`NBodyChargeCloud`](@ref), counting the center as
shell 0.

Successive shells shrink by a factor `e`, so that the outer ones read as a halo around the center
rather than hiding it.
"""
charge_cloud_shell_markersize(markersize::Real, shell::Integer) = markersize * exp(-shell)

"""
    potential_plot_title(name, cross_section, value, units; sigdigits = 3) -> String

The title of a cross-section plot: `"Electric Potential @ φ = 0.0°"`, `"Point Type Map @ z = 40.0 mm"`.

`value` is the cross-section position as `get_crosssection_idx_and_value` returns it — degrees for a
`:φ` cross section, internal length units otherwise — and `units` is the unit it is *shown* in.
`sigdigits` is `3` for the scalar fields and `2` for the electric field magnitude, which is the
(arbitrary) difference the two Plots recipes have always had.

Shared by both backends so the two cannot drift apart; the Makie figure-level wrappers set the same
string on the `Axis` that the Plots recipes set as `title`. They differ in one argument: the Plots
recipes pass the unit the *caller* used (a bare number meaning internal units, so `"@ z = 0.0375 m"`),
while the Makie wrappers pass `display_unit(units)`, so the title reads `"@ z = 37.5 mm"` and agrees
with the axis labels. `test/parity/gate.jl` converts the golden title before comparing.
"""
function potential_plot_title(name::AbstractString, cross_section::Symbol, value::Real,
                              units::Unitful.Units; sigdigits::Int = 3)
    q = value * (cross_section == :φ ? u"°" : internal_length_unit)
    "$name @ $(cross_section) = $(round(units, Float64(uconvert(units, q)), sigdigits=sigdigits))"
end


_get_potential_plot_information(sp::ScalarPotential, ::Bool) = _get_potential_plot_information(sp) # ignore point type bool for all other potentials

function _get_potential_plot_information(p::ElectricPotential)::Tuple{Symbol, Tuple{Quantity, Quantity}, String, Unitful.Units}
    return (:viridis, extrema(p.data).*internal_voltage_unit, "Electric Potential", internal_voltage_unit)
end

function _get_potential_plot_information(::WeightingPotential{T})::Tuple{Symbol, Tuple{T,T}, String, Unitful.Units} where {T}
    return (:viridis, (zero(T), one(T)), "Weighting Potential", Unitful.NoUnits)
end

function _get_potential_plot_information(ecd::EffectiveChargeDensity)::Tuple{Symbol, Tuple{Quantity, Quantity}, String, Unitful.Units}
    return (:inferno, extrema(ecd.data).*u"V/m^2", "Effective Charge Density", u"V/m^2")
end

function _get_potential_plot_information(::PointTypes{T}, all_point_types::Bool)::Tuple{Symbol, Tuple{T,T}, String, Unitful.Units} where {T}
    return (:viridis, (zero(0), (all_point_types ? 31 : 7)*one(T)), "Point Type Map", Unitful.NoUnits)
end

function _get_potential_plot_information(::ImpurityScale{T})::Tuple{Symbol, Tuple{T,T}, String, Unitful.Units} where {T}
    return (:viridis, (zero(T), one(T)), "Impurity Scale", Unitful.NoUnits)
end

# The Plots.jl recipe for the dielectric distribution predates this helper and spells its colour
# scheme, title and (absent) `clims` out inline, so this method is currently used only by the Makie
# extension — which needs a colour range for the `Colorbar` it has to create explicitly (plan §1.4).
# Plots let GR take the limits from the plotted slice; the whole-array extrema used here are the
# same numbers for every detector where a slice contains both the semiconductor and its surroundings,
# and they keep the colour scale comparable between two cross sections of the same detector.
function _get_potential_plot_information(ϵ::DielectricDistribution{T})::Tuple{Symbol, Tuple{T,T}, String, Unitful.Units} where {T}
    return (:inferno, extrema(ϵ.data), "Dielectric Distribution", Unitful.NoUnits)
end

function get_crosssection_idx_and_value(grid::Grid, ::Any, ::Any, ::Any)
    error(ArgumentError, ": The cross section your are trying to plot is not given in the right format.\nPlease pass a number with matching or no units.")
end


## Cross sections of cylindrical grids

function get_crosssection_idx_and_value(grid::CylindricalGrid{T}, r::Union{Real, LengthQuantity, Missing}, φ::Union{Real, AngleQuantity, Missing}, z::Union{Real, LengthQuantity, Missing})::Tuple{Symbol,Int,T,Unitful.Units} where {T <: SSDFloat}

    cross_section::Symbol, idx::Int, units::Unitful.Units = if ismissing(φ) && ismissing(r) && ismissing(z)
        return get_crosssection_idx_and_value(grid, r, T(0.0), z)
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        φ_rad::T = T(to_internal_units(φ isa Real ? deg2rad(φ) : φ))
        while !(grid.φ.interval.left <= φ_rad <= grid.φ.interval.right) && grid.φ.interval.right != grid.φ.interval.left
            if φ_rad > grid.φ.interval.right
                φ_rad -= width(grid.φ.interval)
            elseif φ_rad < grid.φ.interval.left
                φ_rad += width(grid.φ.interval)
            end
        end
        :φ, searchsortednearest(grid.φ, φ_rad), φ isa Real ? u"°" : unit(φ)
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(grid.r, T(to_internal_units(r))), r isa Real ? internal_length_unit : unit(r)
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(grid.z, T(to_internal_units(z))), z isa Real ? internal_length_unit : unit(z)
    else
        error(ArgumentError, ": Only one of the keywords `r, φ, z` is allowed.")
    end

    value::T = if cross_section == :φ
        rad2deg(grid.φ[idx])
    elseif cross_section == :r
        grid.r[idx]
    elseif cross_section == :z
        grid.z[idx]
    end

    cross_section, idx, value, units
end


## Cross sections of Cartesian grids

function get_crosssection_idx_and_value(grid::CartesianGrid3D{T}, x::Union{Real, LengthQuantity, Missing}, y::Union{Real, LengthQuantity, Missing}, z::Union{Real, LengthQuantity, Missing})::Tuple{Symbol,Int,T,Unitful.Units} where {T <: SSDFloat}

    cross_section::Symbol, idx::Int, units::Unitful.Units = if ismissing(x) && ismissing(y) && ismissing(z)
        return get_crosssection_idx_and_value(grid, T(0.0), y, z)
    elseif !ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(grid.x, T(to_internal_units(x))), x isa Real ? internal_length_unit : unit(x)
    elseif ismissing(x) && !ismissing(y) && ismissing(z)
        :y, searchsortednearest(grid.y, T(to_internal_units(y))), y isa Real ? internal_length_unit : unit(y)
    elseif ismissing(x) && ismissing(y) && !ismissing(z)
        :z, searchsortednearest(grid.z, T(to_internal_units(z))), z isa Real ? internal_length_unit : unit(z)
    else
        error(ArgumentError, ": Only one of the keywords `x, y, z` is allowed.")
    end
    value::T = grid[cross_section][idx]
    cross_section, idx, value, units
end


# The title and colorbar title of a field-magnitude plot. Shared with the Makie extension, and kept
# here rather than in a recipe file so they survive the removal of the Plots.jl recipes.
const EFIELD_PLOT_NAME = "Electric Field (Magn.)"
const EFIELD_COLORBAR_NAME = "Electric Field Strength"


## Electric field lines
#
# Which cross section, where charges are spawned, where they drift, and how a path projects onto the
# plane are all backend-agnostic, so both recipes read them from here. The algorithm is the drift
# integration the Plots recipe has always used; `streamplot` is a phase-2 replacement, not a port.

"""
    fieldline_cross_section(sim; x, y, z, r, φ) -> (dim_symbol, value, units)

The cross section a field-line plot is drawn on.

`value` is in internal units (metres, or radians for `:φ`) and `units` is the unit the caller gave it
in, which is what the title is printed in. Passing none selects the middle `x` slice of a Cartesian
grid and `φ = 0` of a cylindrical one. Throws `ArgumentError` if more than one is given.
"""
function fieldline_cross_section(sim::Simulation; x = missing, y = missing, z = missing,
                                 r = missing, φ = missing)
    S = get_coordinate_system(sim.electric_field.grid)
    T = get_precision_type(sim.detector)
    given = S == Cartesian ? [x, y, z] : [r, φ, z]
    names = S == Cartesian ? [:x, :y, :z] : [:r, :φ, :z]

    if all(ismissing, given)
        S == Cartesian || return (:φ, zero(T), u"°")
        ax = sim.electric_field.grid[1]
        return (:x, T(ax[div(length(ax), 2)]), internal_length_unit)
    end
    count(!ismissing, given) == 1 || throw(ArgumentError(
        "Only one keyword for a certain dimension is allowed. Please choose one of " *
        (S == Cartesian ? "'x', 'y', 'z'." : "'r', 'φ', 'z'.")))

    i = findfirst(!ismissing, given)
    name, v = names[i], given[i]
    units = v isa Real ? (name == :φ ? u"°" : internal_length_unit) : unit(v)
    (name, T(to_internal_units(name == :φ && v isa Real ? deg2rad(v) : v)), units)
end

"""
    fieldline_title(dim_symbol, value, units; full_det = false) -> String

The title of a field-line plot: `"Electric Field Lines @ φ = 0.0°"`.

`value` is in internal units and `units` is the unit it is shown in, as
[`potential_plot_title`](@ref) takes them. `full_det` appends the angle mirrored onto the left half.
"""
function fieldline_title(dim_symbol::Symbol, value::Real, units::Unitful.Units;
                         full_det::Bool = false)
    u = dim_symbol == :φ ? u"rad" : internal_length_unit
    at = round(units, Float64(uconvert(units, value * u)), sigdigits = 3)
    title = "Electric Field Lines @ $(dim_symbol) = $(at)"
    full_det || return title
    mirrored = round(units, Float64(uconvert(units, (value + π) % 2π * u"rad")), sigdigits = 3)
    title * "\n(= $(mirrored) on left side)"
end

# `full_det` mirrors the right half onto the left, which only means anything for a φ cross section of
# a cylindrical grid: no other cut has two halves to mirror.
fieldline_full_det(dim_symbol::Symbol, full_det::Bool) = full_det && dim_symbol == :φ

"""
    fieldline_spawn_positions(sim, dim_symbol, value; sampling, offset, spacing, skip_contact)

Points just inside the semiconductor where a field line starts.

The contact surfaces are intersected with the sample lines of `get_sample_lines`, and each hit is
moved `offset` along the inward surface normal so that the charge is in the semiconductor rather than
in the contact, where it would not drift. `spacing = 2` keeps every second point.

`skip_contact` is the contact *not* sampled, normally the core, so that lines start on one polarity
only.
"""
function fieldline_spawn_positions(sim::Simulation, dim_symbol::Symbol, value::Real;
                                   sampling = 2u"mm", offset = 0.5u"mm", spacing::Int = 1,
                                   skip_contact::Int = 1)
    T = get_precision_type(sim.detector)
    positions = CartesianPoint{T}[]
    lines = get_sample_lines(dim_symbol, T(value), sim.electric_field.grid,
                             T(ustrip(to_internal_units(sampling))))
    Δ = T(ustrip(to_internal_units(offset)))
    # The loop nesting is contact, then sample line, then surface. It fixes the order of the points,
    # which `spacing` then thins by index, so reordering it silently changes which lines are drawn.
    for c in sim.detector.contacts
        c.id == skip_contact && continue
        surfs = ConstructiveSolidGeometry.surfaces(c.geometry)
        for l in lines, surf in surfs
            for pt in ConstructiveSolidGeometry.intersection(surf, l)
                pt in c || continue
                inside = pt - Δ * normalize(ConstructiveSolidGeometry.normal(surf, pt))
                if inside in sim.detector && !(inside in sim.detector.contacts)
                    push!(positions, inside)
                end
            end
        end
    end
    spacing <= 1 && return positions
    positions[findall(i -> ((spacing - 1) + i) % spacing == 0, eachindex(positions))]
end

"""
    fieldline_paths(sim, dim_symbol, value; kwargs...) -> Vector{Vector{CartesianPoint}}

The field lines through the given cross section, as drifted paths.

Each spawn position is drifted both along and against the field, so a line is traced out in both
directions from the contact it starts on. The drift uses `ElectricFieldChargeDriftModel`, i.e. a
velocity parallel to **E**, which is what makes the paths field lines rather than carrier
trajectories.

A path that has not reached a contact after `max_nsteps` stops there; raise it, or lower `Δt`, if the
lines end mid-detector.
"""
function fieldline_paths(sim::Simulation, dim_symbol::Symbol, value::Real;
                         max_nsteps::Int = 5000, sampling = 2u"mm", offset = 0.5u"mm",
                         spacing::Int = 1, Δt = 1u"ns", skip_contact::Int = 1)
    T = get_precision_type(sim.detector)
    spawn = fieldline_spawn_positions(sim, dim_symbol, value;
                                      sampling, offset, spacing, skip_contact)
    out = Vector{Vector{CartesianPoint{T}}}()
    isempty(spawn) && return out
    forward = interpolated_vectorfield(sim.electric_field.data, sim.electric_field.grid)
    backward = interpolated_vectorfield(sim.electric_field.data .* -1, sim.electric_field.grid)
    for field in (forward, backward)
        paths = Array{CartesianPoint{T}, 2}(undef, length(spawn), max_nsteps)
        last_step = _drift_charge!(paths, Vector{T}(undef, max_nsteps), sim.detector,
                                   sim.point_types, sim.electric_potential.grid,
                                   CartesianPoint.(spawn), ones(T, length(spawn)),
                                   T(to_internal_units(Δt)), field, Electron,
                                   geometry_check = false, verbose = false)
        for i in axes(paths, 1)
            push!(out, collect(paths[i, 1:last_step]))
        end
    end
    out
end

"""
    fieldline_projection(path, dim_symbol, value, S) -> (u, v)

A field line projected onto its cross-section plane, as two vectors of bare numbers.

The units depend on the cut: a `:z` cross section of a cylindrical grid is polar and gives
(angle in radians, radius in metres); an `:r` cross section gives (angle in radians, height in
metres); every other cut gives two lengths in metres.
"""
function fieldline_projection(path::AbstractVector{<:CartesianPoint}, dim_symbol::Symbol,
                              value::Real, S::Type{<:AbstractCoordinateSystem})
    if dim_symbol == :φ
        # Project onto the φ plane: the signed distance along the cut direction, and z.
        return ([p[1] * cos(value) + p[2] * sin(value) for p in path], [p[3] for p in path])
    elseif dim_symbol == :x
        return ([p[2] for p in path], [p[3] for p in path])
    elseif dim_symbol == :y
        return ([p[1] for p in path], [p[3] for p in path])
    elseif dim_symbol == :z
        S == Cylindrical || return ([p[1] for p in path], [p[2] for p in path])
        cyl = CylindricalPoint.(path)
        return ([p[2] for p in cyl], [p[1] for p in cyl])
    else # :r — the cylindrical shell at that radius, unrolled into (φ, z)
        cyl = CylindricalPoint.(path)
        return ([p[2] for p in cyl], [p[3] for p in cyl])
    end
end


## Field-line sampling lines (used by electric field line plotting)

function get_sample_lines(dim_symbol::Symbol, v::T, grid::CartesianGrid3D{T}, sampling::T)::Vector{ConstructiveSolidGeometry.Line} where {T}

    xrange = range(round.(endpoints(grid.x.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    yrange = range(round.(endpoints(grid.y.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    zrange = range(round.(endpoints(grid.z.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)

    return if dim_symbol == :x
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(v,y,0), CartesianVector{T}(0,0,1)) for y in yrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(v,0,z), CartesianVector{T}(0,1,0)) for z in zrange]
        )
    elseif dim_symbol == :y
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(x,v,0), CartesianVector{T}(0,0,1)) for x in xrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,v,z), CartesianVector{T}(1,0,0)) for z in zrange]
        )

    else # dim_symbol == :z
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(x,0,v), CartesianVector{T}(0,1,0)) for x in xrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,y,v), CartesianVector{T}(1,0,0)) for y in yrange]
        )
    end
end


function get_sample_lines(dim_symbol::Symbol, v::T, grid::CylindricalGrid{T}, sampling::T)::Vector{ConstructiveSolidGeometry.Line} where {T}

    φsampling = 2π * sampling / (dim_symbol == :r ? v : grid.r.interval.right) # or use sampling together with the unit ?
    rrange = range(round.((-grid.r.interval.right,grid.r.interval.right)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    φrange = range(round.(endpoints(grid.φ.interval)./φsampling, (RoundDown, RoundUp)).*φsampling..., step = φsampling)
    zrange = range(round.(endpoints(grid.z.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)

    return if dim_symbol == :φ
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(r*cos(v),r*sin(v),0), CartesianVector{T}(0,0,1)) for r in rrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,0,z), CartesianVector{T}(cos(v),sin(v),0)) for z in zrange]
        )
    elseif dim_symbol == :z
        [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,0,v), CartesianVector{T}(cos(φ), sin(φ),0)) for φ in φrange]
    else # dim_symbol == :r
        [ConstructiveSolidGeometry.Line(CartesianPoint{T}(v*cos(φ),v*sin(φ),0), CartesianVector{T}(0,0,1)) for φ in φrange]
    end
end
