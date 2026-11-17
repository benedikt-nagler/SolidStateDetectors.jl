# Makie recipe and figure-level wrapper for the initial distribution of charge in an event: a
# `PlatonicSolid` (of which `PointCharge` is the one-vertex case), an `NBodyChargeCloud`, and the
# `AbstractParticleSource` that deposited the energy. 11 of the ~142 parity cases.
#
# A cloud is a set of point charges in groups: the center, then one group per shell, each drawn
# smaller than the last, with the edges of the shell's platonic solid dotted in between. The Plots
# recipes emit one scatter series per group and one line series per edge path, which for the default
# two-dodecahedron cloud is 2 + 2 × 10 series. Here `_cloud_groups` flattens the same structure into a
# single `scatter!` carrying a marker size *per point* and a single `lines!` carrying every edge path
# with NaN vertices between them, so the whole cloud is two plot objects whatever it is made of.
#
# A particle source is not a charge cloud, but it is the same picture from one step earlier - where
# the energy came from rather than where it landed - and it is what the Plots recipes group it with,
# so `chargecloudplot` draws it too.

"""
    chargecloudplot(cloud; kwargs...)
    chargecloudplot!(ax, cloud; kwargs...)

Plot the initial charge distribution of an event: an [`NBodyChargeCloud`](@ref), a single
[`PointCharge`](@ref) or any other `PlatonicSolid` shell, or the `AbstractParticleSource` the energy
came from.

The non-mutating form returns a `FigureAxisPlot` carrying the axis labels and the legend; the
mutating form draws into an existing axis and sets neither. `plot(cloud)` routes here.

Coordinates are drawn in millimetres, so a cloud can be drawn into the same axis as a
[`detectorplot`](@ref) of the detector it sits in.

## Keywords
* `connect`: draw the edges of each shell. Off is the readable choice for a cloud of many shells.
* `color`: colour of the carriers.
* `markersize`: size of the center carrier; each shell is drawn a factor `e` smaller than the one
  inside it.
* `linewidth`, `linestyle`: the shell edges.
* `length`: for a directed source, how far its ray or emission cone reaches, in internal units.

A source emitting isotropically (`opening_angle = 180°`, which is what the two-argument
`MonoenergeticSource` constructor gives) is drawn as its position alone: there is no direction to
draw.
"""
@recipe ChargeCloudPlot (cloud,) begin
    "Draw the edges of each shell."
    connect = true
    "Colour of the carriers; `Makie.automatic` is blue."
    color = Makie.automatic
    "Size of the center carrier; shells shrink by a factor `e` each."
    markersize = Makie.automatic
    "Width of the shell edges."
    linewidth = Makie.automatic
    "Style of the shell edges."
    linestyle = :dot
    "Reach of a directed source's ray or emission cone, in internal units."
    length = 0.01
end

const ChargeCloudLike = Union{AbstractChargeCloud, AbstractParticleSource}

for GT in (:AbstractChargeCloud, :AbstractParticleSource)
    @eval Makie.plottype(::$GT) = ChargeCloudPlot
    @eval Makie.convert_arguments(::Type{<:ChargeCloudPlot}, c::$GT) = (c,)
end

Makie.preferred_axis_type(::ChargeCloudPlot) = Makie.Axis3


## ----------------------------------------------------------------------------------------------
## Decomposition
## ----------------------------------------------------------------------------------------------

"""
    _cloud_groups(cloud, markersize) -> Vector{NamedTuple}

The carriers of a charge cloud, grouped by the marker size they are drawn at.

Each group carries its `points`, its `markersize`, and the `connections` between its points as index
paths into the group (see `charge_cloud_connections`). A `PlatonicSolid` is one group; an
`NBodyChargeCloud` is its center plus one group per shell.
"""
_cloud_groups(s::PlatonicSolid, markersize::Real) =
    [(points = collect(s.locations), markersize = Float64(markersize),
      connections = charge_cloud_connections(s))]

function _cloud_groups(nbc::NBodyChargeCloud{T, N}, markersize::Real) where {T, N}
    pts = nbc.locations
    groups = [(points = [pts[1]], markersize = Float64(markersize), connections = Vector{Int}[])]
    first_vertex = 1
    for shell in Base.OneTo(N)
        n = _shell_size(nbc.shell_structure, shell)
        push!(groups, (points = pts[(first_vertex + 1):(first_vertex + n)],
                       markersize = charge_cloud_shell_markersize(markersize, shell),
                       connections = _shell_connections(nbc.shell_structure)))
        first_vertex += n
    end
    groups
end

# The two `NBodyChargeCloud` constructors give `shell_structure` different meanings: the shell shape
# itself, or - for the `N`-charges-per-shell constructor, which distributes points over a sphere with
# no structure to connect - the number of charges in each shell.
_shell_size(SH::Type{<:PlatonicSolid}, ::Int) = get_vertices(SH)
_shell_size(sizes::NTuple{<:Any, Int}, shell::Int) = sizes[shell]

_shell_connections(SH::Type{<:PlatonicSolid}) = charge_cloud_connections(SH)
_shell_connections(::NTuple{<:Any, Int}) = Vector{Int}[]


## ----------------------------------------------------------------------------------------------
## The recipe
## ----------------------------------------------------------------------------------------------

function Makie.plot!(p::ChargeCloudPlot{<:Tuple{<:AbstractChargeCloud}})
    cloud = p[1][]
    color = p.color[] === Makie.automatic ? SSD_COLORS.charge_cloud : p.color[]
    ms = p.markersize[] === Makie.automatic ? _default_markersize(cloud) : p.markersize[]
    groups = _cloud_groups(cloud, ms)

    pts, sizes = Makie.Point3f[], Float64[]
    for g in groups
        append!(pts, point_coordinates(g.points))
        append!(sizes, fill(g.markersize, length(g.points)))
    end
    Makie.scatter!(p, pts; color, markersize = sizes, strokewidth = 0)

    if p.connect[]
        edges = _cloud_edges(groups)
        isempty(edges) || Makie.lines!(p, edges; color, linestyle = p.linestyle[],
                                       linewidth = p.linewidth[] === Makie.automatic ?
                                           SSD_LINEWIDTHS.charge_cloud : p.linewidth[])
    end
    pad_flat_extent!(p)
    return p
end

# A bare shell is drawn at the size the Plots `PlatonicSolid` recipe gives it; a whole cloud at the
# larger size its `NBodyChargeCloud` recipe gives the center, since its shells shrink from there.
_default_markersize(::AbstractChargeCloud) = SSD_MARKERSIZES.charge_cloud
_default_markersize(::NBodyChargeCloud) = SSD_MARKERSIZES.charge_cloud_center

# Every edge path of every group, end to end, with a non-finite vertex between paths so that Makie
# lifts the pen rather than joining one path to the next.
function _cloud_edges(groups)
    pts = Makie.Point3f[]
    for g in groups
        vertices = point_coordinates(g.points)
        for path in g.connections
            isempty(pts) || push!(pts, Makie.Point3f(NaN32, NaN32, NaN32))
            append!(pts, vertices[path])
        end
    end
    pts
end


## ----------------------------------------------------------------------------------------------
## Particle sources
## ----------------------------------------------------------------------------------------------

function Makie.plot!(p::ChargeCloudPlot{<:Tuple{<:AbstractParticleSource}})
    source = p[1][]
    _plot_emission!(p, source)
    Makie.scatter!(p, point_coordinates(source.position);
                   color = p.color[] === Makie.automatic ? SSD_COLORS.particle_source : p.color[],
                   markersize = p.markersize[] === Makie.automatic ?
                       SSD_MARKERSIZES.particle_source : p.markersize[],
                   strokecolor = :black, strokewidth = 1)
    pad_flat_extent!(p)
    return p
end

# What the source emits into: a ray for a perfectly directed source, a cone for one with an opening
# angle up to 90°, and nothing at all beyond that - an isotropic source has no direction to draw, and
# a cone wider than a hemisphere would hide the geometry around it rather than describe the source.
function _plot_emission!(p, source)
    hasproperty(source, :opening_angle) || return nothing
    reach = p.length[]
    if iszero(source.opening_angle)
        tip = source.position + reach * normalize(source.direction)
        Makie.lines!(p, point_coordinates([source.position, tip]);
                     color = SSD_COLORS.emission, linewidth = 2)
    elseif source.opening_angle <= 90u"°"
        csgplot!(p, _emission_cone(source, reach);
                 color = SSD_COLORS.emission, linewidth = 0, pad_flat = false)
    end
    return nothing
end

# The cone the source emits into, with its apex at the source and its axis along `direction`. Built
# from an orthonormal frame around `direction`, since a `Cone` is defined along its own z axis.
function _emission_cone(source, reach)
    d = normalize(source.direction)
    a = normalize(d × (abs(d.x) == 1 ? CartesianVector(0, 1, 0) : CartesianVector(1, 0, 0)))
    b = normalize(a × d)
    rot = hcat(a, b, d)
    half_height = reach * cos(source.opening_angle) / 2
    Cone(r = ((0, 0), (0, reach * sin(source.opening_angle))), hZ = half_height,
         origin = source.position + rot * CartesianVector(0, 0, half_height),
         rotation = rot)
end


## ----------------------------------------------------------------------------------------------
## Figure-level wrapper
## ----------------------------------------------------------------------------------------------

# The legend entry Plots produced from the label of the first series, as an explicit element (§1.3):
# one swatch naming the object, from `plot_label`, which the Plots recipes call too.
function _cloud_legend_entries(cloud; color = Makie.automatic, kwargs...)
    swatch = color === Makie.automatic ? _default_color(cloud) : color
    (Makie.LegendElement[charge_cloud_legend_element(swatch)], String[plot_label(cloud)])
end

_default_color(::AbstractChargeCloud) = SSD_COLORS.charge_cloud
_default_color(::AbstractParticleSource) = SSD_COLORS.particle_source

function chargecloudplot(gp::Makie.GridPosition, cloud::ChargeCloudLike;
                         axis = (;), legend = (;), kwargs...)
    gl = panel_layout(gp)
    ax = ssd_axis(gl[1, 1]; geometry_axis_setup().axis..., axis...)
    plt = chargecloudplot!(ax, cloud; kwargs...)
    finish_geometry_axis!(ax, plt, axis)
    legend === false || add_legend!(gl, _cloud_legend_entries(cloud; kwargs...)...;
                                    corner = legend_corner(ax, plt), legend...)
    Makie.AxisPlot(ax, plt)
end

function chargecloudplot(cloud::ChargeCloudLike; figure = (;), kwargs...)
    fig = Makie.Figure(; figure...)
    ap = chargecloudplot(fig[1, 1], cloud; kwargs...)
    Makie.FigureAxisPlot(fig, ap.axis, ap.plot)
end

Makie.plot(cloud::ChargeCloudLike; kwargs...) = chargecloudplot(cloud; kwargs...)
