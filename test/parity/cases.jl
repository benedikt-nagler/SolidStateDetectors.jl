# The plot-parity case matrix.
#
# Every case is one row of the matrix: it is extracted into a golden spec from Plots.jl
# (`devscripts/extract_plots_spec.jl`), asserted against Makie in `test/test_plotting.jl`, and
# rendered as a side-by-side pair by `devscripts/plot_gallery.jl`. The matrix is generated
# programmatically from three fixtures — a cylindrical simulation, a Cartesian one, and one drifted
# event — so it stays a data structure rather than 130 hand-written blocks.
#
# This file must load without either plotting backend: it only builds objects and records *how* they
# are to be plotted.

module PlotParityCases

using SolidStateDetectors
using SolidStateDetectors: ConstructiveSolidGeometry as CSG
using SolidStateDetectors: Simulation, Event, SSD_examples, simulate!
using Unitful

export Case, cases, fixtures, case_ids, group_ids, call_string, NOTE_R_ALIAS

const T = Float32

"""
    Case

One row of the parity matrix.

* `id`: filename-safe unique identifier; also the name of the golden JSON.
* `group`: coarse grouping, matching the plan's §3.4 table. Used to filter runs.
* `kind`: how the case is called.
  - `:plot` — `plot(args...; kwargs...)`, the overwhelming majority.
  - `:fieldlines` — the `plot_electric_fieldlines(sim; kwargs...)` user plot.
  - `:ehcontribution` — the `plot_electron_and_hole_contribution(evt, sim, id; kwargs...)` user plot.
* `args`: a thunk returning the positional arguments. A thunk, not the objects themselves, so that
  building the matrix is cheap and a case that cannot be constructed fails only that case.
* `kwargs`: keyword arguments, as passed to Plots.
* `slow`: cases that drift charges and take seconds to minutes each.
* `note`: non-empty when the Plots output of this case is *not* a valid parity reference. Filled in
  automatically for the cases hit by `NOTE_R_ALIAS`.

`kwargs` are the *Plots* keywords. Most are semantic and carry over unchanged (`r`, `φ`, `z`,
`full_det`, `contours_equal_potential`, …); `seriestype` is the exception — for detector and CSG
cases it selects the representation (`:csg`, `:slice`, `:wireframe`, `:mesh3d`, `:samplesurface`),
and the Makie side maps it onto its own entry point when those steps land.
"""
struct Case
    id::String
    group::String
    kind::Symbol
    args::Function
    kwargs::NamedTuple
    slow::Bool
    note::String
end

Case(id, group, args, kwargs = NamedTuple(); kind = :plot, slow = false, note = "") =
    Case(String(id), String(group), kind, args, kwargs, slow, note)

"""
    call_string(c::Case, args = nothing)

Human-readable rendering of the call, stored in the spec and shown in the gallery. Pass the already
built `args` to have the argument types spelled out; `call_string` never builds them itself.
"""
function call_string(c::Case, args = nothing)
    fname = c.kind === :plot ? "plot" :
            c.kind === :fieldlines ? "plot_electric_fieldlines" :
            "plot_electron_and_hole_contribution"
    argstr = args === nothing ? "…" : join(string.(nameof.(typeof.(args))), ", ")
    kwstr = join(["$k = $(repr(v))" for (k, v) in pairs(c.kwargs)], ", ")
    "$fname($argstr$(isempty(kwstr) ? "" : "; " * kwstr))"
end

## ----------------------------------------------------------------------------------------------
## Fixtures
## ----------------------------------------------------------------------------------------------

const _SIMULATION_FIXTURES = Ref{Any}(nothing)
const _EVENT_FIXTURE = Ref{Any}(nothing)

"""
    fixtures(; refresh = false, with_event = false)

The simulated objects the matrix uses, computed once per session:

* `cyl`: `:InvertedCoax`, a cylindrical-grid simulation.
* `car`: `:CGD`, a Cartesian-grid simulation.
* `evt`: one event drifted in `cyl`, with waveforms. It is built only when `with_event = true`.

`refinement_limits = missing` keeps the grids at their unrefined default: fast, and — more
importantly — deterministic, so the golden specs do not move when the refinement heuristics change.
"""
function fixtures(; refresh::Bool = false, with_event::Bool = false)
    if refresh
        _SIMULATION_FIXTURES[] = nothing
        _EVENT_FIXTURE[] = nothing
    end
    if _SIMULATION_FIXTURES[] === nothing
        cyl = Simulation{T}(SSD_examples[:InvertedCoax])
        simulate!(cyl, device_array_type = Array, refinement_limits = missing, verbose = false)

        car = Simulation{T}(SSD_examples[:CGD])
        simulate!(car, device_array_type = Array, refinement_limits = missing, verbose = false)
        _SIMULATION_FIXTURES[] = (; cyl, car)
    end
    sims = _SIMULATION_FIXTURES[]
    if with_event && _EVENT_FIXTURE[] === nothing
        evt = Event([CartesianPoint{T}(0.02, 0.01, 0.05)])
        simulate!(evt, sims.cyl)
        _EVENT_FIXTURE[] = evt
    end
    (; sims..., evt = with_event ? _EVENT_FIXTURE[] : missing)
end

## ----------------------------------------------------------------------------------------------
## The matrix
## ----------------------------------------------------------------------------------------------

# `r` is a Plots *attribute alias* for `:rotation` (`Plots._keyAliases[:r] === :rotation`), so
# `plot(potential; r = 0.02)` is rewritten before the recipe's keyword extraction runs and the recipe
# never sees it — it silently falls back to its no-selector default, the `φ = 0` cross-section. Every
# `r =` cross-section is therefore unreachable through Plots today, and the specs extracted for those
# cases describe the fallback, not the intended plot. They are kept in the matrix because Makie has
# no such alias and must implement the cross-section properly; they just carry no usable reference.
const NOTE_R_ALIAS = "No valid Plots reference: `r` is a Plots alias for `:rotation`, so the \
cross-section keyword never reaches the recipe and the recorded plot is the `φ = 0` fallback. \
Verify this case at tier 2 (visual review) instead of against the golden spec."

# `plot(::CartesianVector)` throws under Plots 1.41: the `Val{:vector}` series-type recipe in
# `src/ConstructiveSolidGeometry/plotting/PointsAndVectors/Vectors.jl` runs only *after* Plots has
# converted the positional arguments into series data, and `(CartesianPoint, CartesianVector)` has no
# such conversion. So the hand-built 13-point arrowhead is unreachable and there is no reference
# image to match. Makie draws vectors with `arrows3d!` (plan §2), so the case stays in the matrix.
const NOTE_VECTOR_BROKEN = "No valid Plots reference: `plot(::CartesianVector)` errors under \
Plots 1.41 (the `:vector` series-type recipe is never reached). Verify at tier 2 instead."

# Cross-section keywords per coordinate system, paired with the ASCII name they get in a case id
# (case ids become filenames, so `φ` is spelled `phi`). The values are grid ticks of the fixtures;
# the recipes snap to the nearest tick anyway, so the exact values are not critical — but they are
# fixed here so the golden specs are reproducible.
const CYL_SLICES = (("r", (r = 0.02,)), ("phi", (φ = 0,)), ("z", (z = 0.04,)))
const CAR_SLICES = (("x", (x = 0,)), ("y", (y = 0,)), ("z", (z = 0,)))

# Scalar fields of a Simulation and the name they get in a case id.
const POTENTIAL_FIELDS = (
    electricpotential   = sim -> sim.electric_potential,
    weightingpotential  = sim -> sim.weighting_potentials[1],
    chargedensity       = sim -> sim.q_eff_imp,
    pointtypes          = sim -> sim.point_types,
    impurityscale       = sim -> sim.imp_scale,
    dielectric          = sim -> sim.ϵ_r,
)

function _potential_cases()
    cs = Case[]
    for (csname, slices, simfield) in (("cyl", CYL_SLICES, :cyl), ("car", CAR_SLICES, :car))
        for (fieldname, getter) in pairs(POTENTIAL_FIELDS)
            for (slicename, slicekw) in slices
                push!(cs, Case("potentials_$(csname)_$(fieldname)_$(slicename)", "potentials",
                               () -> (getter(getfield(fixtures(), simfield)),), slicekw))
            end
        end
    end

    # Variants. The dielectric distribution takes no extra keywords, so it is left out of these.
    cyl = () -> fixtures().cyl
    car = () -> fixtures().car
    push!(cs,
        Case("potentials_cyl_electricpotential_phi_contours", "potentials",
             () -> (cyl().electric_potential,), (φ = 0, contours_equal_potential = true)),
        Case("potentials_cyl_electricpotential_phi_fulldet", "potentials",
             () -> (cyl().electric_potential,), (φ = 0, full_det = true)),
        Case("potentials_cyl_electricpotential_phi_fulldet_contours", "potentials",
             () -> (cyl().electric_potential,), (φ = 0, full_det = true, contours_equal_potential = true)),
        Case("potentials_cyl_weightingpotential_phi_contours", "potentials",
             () -> (cyl().weighting_potentials[1],), (φ = 0, contours_equal_potential = true)),
        Case("potentials_cyl_weightingpotential_phi_fulldet", "potentials",
             () -> (cyl().weighting_potentials[1],), (φ = 0, full_det = true)),
        Case("potentials_cyl_pointtypes_phi_all", "potentials",
             () -> (cyl().point_types,), (φ = 0, all_point_types = true)),
        Case("potentials_cyl_pointtypes_r_all", "potentials",
             () -> (cyl().point_types,), (r = 0.02, all_point_types = true)),
        Case("potentials_cyl_chargedensity_phi_fulldet", "potentials",
             () -> (cyl().q_eff_imp,), (φ = 0, full_det = true)),
        Case("potentials_car_electricpotential_y_contours", "potentials",
             () -> (car().electric_potential,), (y = 0, contours_equal_potential = true)),
        Case("potentials_car_weightingpotential_y_contours", "potentials",
             () -> (car().weighting_potentials[1],), (y = 0, contours_equal_potential = true)),
        Case("potentials_car_pointtypes_z_all", "potentials",
             () -> (car().point_types,), (z = 0, all_point_types = true)),
        # A φ cross-section given with an explicit unit exercises the other branch of
        # `get_crosssection_idx_and_value` and shows up in the title and the tick labels.
        Case("potentials_cyl_electricpotential_phi_units", "potentials",
             () -> (cyl().electric_potential,), (φ = 30u"°",)),
        Case("potentials_cyl_electricpotential_r_units", "potentials",
             () -> (cyl().electric_potential,), (r = 20u"mm",)),
    )
    cs
end

function _efield_cases()
    cs = Case[]
    cyl = () -> fixtures().cyl
    car = () -> fixtures().car
    for (slicename, slicekw) in CYL_SLICES
        push!(cs, Case("efield_cyl_$(slicename)", "efield", () -> (cyl().electric_field,), slicekw))
    end
    for (slicename, slicekw) in CAR_SLICES
        push!(cs, Case("efield_car_$(slicename)", "efield", () -> (car().electric_field,), slicekw))
    end
    push!(cs,
        Case("efield_cyl_phi_contours", "efield",
             () -> (cyl().electric_field,), (φ = 0, contours_equal_potential = true)),
        Case("efield_cyl_phi_fulldet", "efield",
             () -> (cyl().electric_field,), (φ = 0, full_det = true)),
        Case("efield_cyl_phi_fulldet_contours", "efield",
             () -> (cyl().electric_field,), (φ = 0, full_det = true, contours_equal_potential = true)),
        Case("efield_cyl_r_contours", "efield",
             () -> (cyl().electric_field,), (r = 0.02, contours_equal_potential = true)),
        # The Cartesian electric-field recipe has no `full_det` keyword — contours only.
        Case("efield_car_y_contours", "efield",
             () -> (car().electric_field,), (y = 0, contours_equal_potential = true)),
        Case("efield_car_z_contours", "efield",
             () -> (car().electric_field,), (z = 0, contours_equal_potential = true)),
    )
    cs
end

# Field lines drift charges through the field for every spawn position, so each of these takes
# seconds to minutes; `sampling` is kept coarse to hold that down. It has to scale with the detector,
# though: the InvertedCoax is 40 mm in radius and 90 mm tall, the CGD is a 10 mm cube, and a 4 mm
# sampling on the latter spawns no charges at all. The cylindrical sampling is coarser than it needs
# to be for a pretty picture: at 4 mm each case emits ~200 line series and the golden JSON grows past
# a quarter of a megabyte, which is not worth committing for a denser version of the same plot.
#
# `fieldlines_car_x` and `fieldlines_car_default` legitimately draw *nothing*: the CGD contacts are
# slabs perpendicular to x, so no contact surface intersects an x cross-section and no charge is ever
# spawned. They are kept because "draws nothing without erroring" is itself a parity property.
function _fieldline_cases()
    cyl = () -> (fixtures().cyl,)
    car = () -> (fixtures().car,)
    fl(id, argf, kw) = Case(id, "fieldlines", argf, kw; kind = :fieldlines, slow = true)
    [
        fl("fieldlines_cyl_phi",          cyl, (φ = 0, sampling = 8u"mm")),
        fl("fieldlines_cyl_phi_fulldet",  cyl, (φ = 0, sampling = 8u"mm", full_det = true)),
        fl("fieldlines_cyl_phi_spacing",  cyl, (φ = 0, sampling = 8u"mm", spacing = 2)),
        fl("fieldlines_cyl_r",            cyl, (r = 0.02, sampling = 8u"mm")),
        fl("fieldlines_cyl_z",            cyl, (z = 0.04, sampling = 8u"mm")),
        fl("fieldlines_cyl_phi_skip2",    cyl, (φ = 0, sampling = 8u"mm", skip_contact = 2)),
        fl("fieldlines_car_x",            car, (x = 0, sampling = 1u"mm")),
        fl("fieldlines_car_y",            car, (y = 0, sampling = 1u"mm")),
        fl("fieldlines_car_z",            car, (z = 0, sampling = 1u"mm")),
        fl("fieldlines_car_default",      car, (sampling = 1u"mm",)),
        fl("fieldlines_cyl_default",      cyl, (sampling = 8u"mm",)),
        fl("fieldlines_car_y_spacing",    car, (y = 0, sampling = 1u"mm", spacing = 3)),
    ]
end

const DETECTOR_VIEWS = (:csg, :wireframe, :mesh3d, :samplesurface)

function _detector_cases()
    cs = Case[]
    det = () -> (fixtures().cyl.detector,)
    detc = () -> (fixtures().car.detector,)
    for st in DETECTOR_VIEWS
        push!(cs, Case("detector_cyl_$(st)", "detector", det, (seriestype = st,)))
    end
    push!(cs,
        Case("detector_cyl_csg_semiconductor", "detector", det,
             (seriestype = :csg, show_semiconductor = true)),
        Case("detector_cyl_csg_nopassives", "detector", det,
             (seriestype = :csg, show_passives = false)),
        Case("detector_cyl_slice_phi", "detector", det, (seriestype = :slice, φ = 0)),
        Case("detector_cyl_slice_x", "detector", det, (seriestype = :slice, x = 0)),
        Case("detector_cyl_slice_z", "detector", det, (seriestype = :slice, z = 0.04)),
        Case("detector_car_csg", "detector", detc, (seriestype = :csg,)),
        Case("detector_car_wireframe", "detector", detc, (seriestype = :wireframe,)),
        Case("detector_car_slice_y", "detector", detc, (seriestype = :slice, y = 0)),
        # Individual parts: a Contact, the Semiconductor, a Passive.
        Case("detector_cyl_contact1", "detector", () -> (fixtures().cyl.detector.contacts[1],), NamedTuple()),
        Case("detector_cyl_contact2", "detector", () -> (fixtures().cyl.detector.contacts[2],), NamedTuple()),
        Case("detector_cyl_contact1_wireframe", "detector",
             () -> (fixtures().cyl.detector.contacts[1],), (seriestype = :wireframe,)),
        Case("detector_cyl_semiconductor", "detector", () -> (fixtures().cyl.detector.semiconductor,), NamedTuple()),
        Case("detector_cyl_semiconductor_slice_phi", "detector",
             () -> (fixtures().cyl.detector.semiconductor,), (seriestype = :slice, φ = 0)),
        Case("detector_car_contact1", "detector", () -> (fixtures().car.detector.contacts[1],), NamedTuple()),
    )
    cs
end

# CSG primitives are constructed here rather than taken from a detector, so the matrix covers
# primitives that no bundled config uses (torus, ellipsoid, prisms) and stays stable if the example
# configs change.
function _csg_cases()
    cs = Case[]
    prim(id, f, kw = NamedTuple()) = Case("csg_$(id)", "csg", f, kw)

    box    = () -> CSG.Box{T}(hX = 1, hY = 2, hZ = 1)
    cone   = () -> CSG.Cone{T}(r = 1, hZ = 1)
    tube   = () -> CSG.Cone{T}(r = ((T(0.5),), (T(1),)), hZ = 1)
    pcone  = () -> CSG.Cone{T}(r = 1, φ = T(π), hZ = 1)
    torus  = () -> CSG.Torus{T}(r_torus = 2, r_tube = 0.5)
    ellips = () -> CSG.Ellipsoid{T}(r = 1)
    prism  = () -> CSG.RegularPrism{6, T}(r = 1, hZ = 1)

    for (name, f) in (("box", box), ("cone", cone), ("tube", tube), ("partialcone", pcone),
                      ("torus", torus), ("ellipsoid", ellips), ("hexprism", prism))
        push!(cs, prim(name, () -> (f(),)))
        push!(cs, prim("$(name)_wireframe", () -> (f(),), (seriestype = :wireframe,)))
    end
    push!(cs,
        prim("box_mesh3d",       () -> (box(),),  (seriestype = :mesh3d,)),
        prim("box_samplesurface", () -> (box(),), (seriestype = :samplesurface,)),
        prim("box_slice_z",      () -> (box(),),  (seriestype = :slice, z = 0)),
        prim("cone_slice_phi",   () -> (cone(),), (seriestype = :slice, φ = 0)),
        # Boolean combinations.
        prim("union",        () -> (box() + cone(),)),
        prim("intersection", () -> (box() & cone(),)),
        prim("difference",   () -> (box() - cone(),)),
        prim("difference_slice_y", () -> (box() - cone(),), (seriestype = :slice, y = 0)),
        # Surface, line and point primitives.
        prim("conemantle",       () -> (CSG.ConeMantle{T}(r = T(1), hZ = T(1)),)),
        prim("torusmantle",      () -> (CSG.TorusMantle{T}(r_torus = T(2), r_tube = T(0.5)),)),
        prim("ellipsoidmantle",  () -> (CSG.EllipsoidMantle{T}(r = T(1)),)),
        prim("ellipticalsurface", () -> (CSG.EllipticalSurface{T}(r = T(1)),)),
        prim("ellipse",  () -> (CSG.Ellipse{T}(r = T(1)),)),
        prim("point",    () -> (CartesianPoint{T}(1, 2, 3),)),
        prim("points",   () -> ([CartesianPoint{T}(cos(φ), sin(φ), 0) for φ in range(0, 2π, length = 17)],)),
        prim("cylpoint", () -> (CylindricalPoint{T}(1, T(π / 4), 2),)),
        # `plot(::CartesianVector)` errors under Plots 1.41 (see NOTE_VECTOR_BROKEN); the case is
        # kept so the Makie port, which draws it with `arrows3d!`, is covered.
        Case("csg_vector", "csg", () -> (CartesianVector{T}(1, 1, 1),), NamedTuple();
             note = NOTE_VECTOR_BROKEN),
    )
    cs
end

function _driftpath_cases()
    evt = () -> fixtures(with_event = true).evt
    [
        Case("driftpaths_single", "driftpaths", () -> (evt().drift_paths[1],)),
        Case("driftpaths_single_nolabel", "driftpaths", () -> (evt().drift_paths[1],), (showlabel = false,)),
        Case("driftpaths_all", "driftpaths", () -> (evt().drift_paths,)),
    ]
end

function _chargecloud_cases()
    center = CartesianPoint{T}(0, 0, 0)
    shells = (:Tetrahedron, :Octahedron, :Hexahedron, :Icosahedron, :Dodecahedron)
    cs = Case[
        Case("chargecloud_pointcharge", "chargeclouds",
             () -> (SolidStateDetectors.PointCharge(center, T(1)),)),
    ]
    for sh in shells
        shape = getfield(SolidStateDetectors, sh)
        push!(cs, Case("chargecloud_$(lowercase(String(sh)))", "chargeclouds",
                       () -> (NBodyChargeCloud(center, T(1) * u"eV";
                                               number_of_shells = 1, shell_structure = shape),)))
    end
    push!(cs,
        Case("chargecloud_icosahedron_noconnect", "chargeclouds",
             () -> (NBodyChargeCloud(center, T(1) * u"eV"; number_of_shells = 1,
                                     shell_structure = SolidStateDetectors.Icosahedron),),
             (connect = false,)),
        Case("chargecloud_nbody_shells", "chargeclouds",
             () -> (NBodyChargeCloud(center, T(1) * u"eV"; number_of_shells = 2,
                                     shell_structure = SolidStateDetectors.Icosahedron),)),
        # The `N`-charges-per-shell constructor takes a different code path (no shell structure,
        # hence no connecting lines).
        Case("chargecloud_nbody_n", "chargeclouds",
             () -> (NBodyChargeCloud(center, T(1) * u"eV", 40; number_of_shells = 2),)),
        Case("chargecloud_monoenergeticsource", "chargeclouds",
             () -> (SolidStateDetectors.MonoenergeticSource(
                        "gamma", 1u"MeV", CartesianPoint(0.04, 0.0, 0.05)),)),
    )
    cs
end

function _waveform_cases()
    evt = () -> fixtures(with_event = true).evt
    [
        Case("waveform_single", "waveforms", () -> (evt().waveforms[1],)),
        Case("waveform_all", "waveforms", () -> (evt().waveforms,)),
        Case("waveform_ehcontribution", "waveforms",
             () -> (evt(), fixtures().cyl, 1); kind = :ehcontribution),
        Case("waveform_ehcontribution_nsamples", "waveforms",
             () -> (evt(), fixtures().cyl, 1), (n_samples = 300,); kind = :ehcontribution),
        Case("waveform_ehcontribution_contact2", "waveforms",
             () -> (evt(), fixtures().cyl, 2); kind = :ehcontribution),
    ]
end

function _grid_cases()
    [
        Case("grid_cyl", "grids", () -> (fixtures().cyl.electric_potential.grid,)),
        Case("grid_car", "grids", () -> (fixtures().car.electric_potential.grid,)),
    ]
end

_with_notes(cs) = [haskey(c.kwargs, :r) && isempty(c.note) ?
                       Case(c.id, c.group, c.args, c.kwargs; c.kind, c.slow, note = NOTE_R_ALIAS) : c
                   for c in cs]

"""
    cases(; groups = nothing, include_slow = true) -> Vector{Case}

The full matrix, optionally filtered to `groups` (a collection of group-name strings) and with the
charge-drifting field-line cases excluded when `include_slow = false`.
"""
function cases(; groups = nothing, include_slow::Bool = true)
    cs = _with_notes(vcat(_potential_cases(), _efield_cases(), _fieldline_cases(),
                          _detector_cases(), _csg_cases(), _driftpath_cases(),
                          _chargecloud_cases(), _waveform_cases(), _grid_cases()))
    groups === nothing || (cs = filter(c -> c.group in groups, cs))
    include_slow || (cs = filter(c -> !c.slow, cs))
    ids = getfield.(cs, :id)
    length(unique(ids)) == length(ids) ||
        error("duplicate case ids: $(join(unique([i for i in ids if count(==(i), ids) > 1]), ", "))")
    cs
end

case_ids(; kwargs...) = getfield.(cases(; kwargs...), :id)
group_ids() = unique(getfield.(cases(), :group))

end # module
