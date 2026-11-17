# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

module SolidStateDetectorsMakieExt

using Makie

using LinearAlgebra
using Rotations
using Unitful

using SolidStateDetectors
using SolidStateDetectors.ConstructiveSolidGeometry

# The plotting entry points are *declared* in the main module (see the "Makie plotting entry points"
# block in `src/SolidStateDetectors.jl`); `import` — not `using` — is what lets `@recipe` attach its
# methods to those functions instead of creating new ones in this module.
import SolidStateDetectors: electricfieldplot, electricfieldplot!,
    scalarpotentialplot, scalarpotentialplot!,
    csgplot, csgplot!, detectorplot, detectorplot!,
    driftpathplot, driftpathplot!, chargecloudplot, chargecloudplot!,
    waveformplot, waveformplot!, ehcontributionplot, ehcontributionplot!,
    gridplot, gridplot!, fieldlinesplot, fieldlinesplot!

# Non-exported helpers / types reused from the main package.
import SolidStateDetectors: SSDFloat, ElectricField, SolidStateDetector,
    Contact, Semiconductor, Passive, Event, Simulation, Grid,
    ScalarPotential, PointTypes, DielectricDistribution,
    EHDriftPath, AbstractChargeCloud, NBodyChargeCloud, PlatonicSolid, AbstractParticleSource,
    get_crosssection_idx_and_value, get_extended_ticks, get_2π_potential,
    _get_potential_plot_information, potential_plot_title,
    plot_label, electron_label, hole_label,
    electron_contribution_label, hole_contribution_label, sum_contribution_label,
    grid_panel_title, waveform_signal_label,
    get_electron_and_hole_contribution, add_baseline_and_extend_tail,
    fieldline_cross_section, fieldline_title, fieldline_full_det, fieldline_paths,
    fieldline_projection, get_coordinate_system, get_precision_type,
    charge_cloud_connections, charge_cloud_shell_markersize, get_vertices,
    EFIELD_PLOT_NAME, EFIELD_COLORBAR_NAME,
    internal_efield_unit, internal_length_unit, internal_angle_unit, internal_voltage_unit
import SolidStateDetectors.ConstructiveSolidGeometry: Cartesian, Cylindrical,
    AbstractGeometry, AbstractPrimitive, AbstractVolumePrimitive, AbstractSurfacePrimitive,
    AbstractPlanarSurfacePrimitive, AbstractLinePrimitive, AbstractConstructiveGeometry,
    AbstractCoordinatePoint, Mesh, Edge, Ellipse, Cone,
    mesh, surfaces, vertices, edges, sample, primitives, isClosedPrimitive, extremum, get_scale,
    get_label_name, get_crosssection_and_value, is_coplanar, rotate, transform
import StatsBase
import StatsBase: midpoints
import SolidStateDetectors.RadiationDetectorSignals: RDWaveform

# Shared infrastructure, used by every recipe and figure-level wrapper below.
include(joinpath(@__DIR__, "Makie", "theme.jl"))
include(joinpath(@__DIR__, "Makie", "units.jl"))
include(joinpath(@__DIR__, "Makie", "axes.jl"))
include(joinpath(@__DIR__, "Makie", "legends.jl"))
include(joinpath(@__DIR__, "Makie", "meshing.jl"))

include(joinpath(@__DIR__, "Makie", "Potentials.jl"))
include(joinpath(@__DIR__, "Makie", "ElectricField.jl"))
include(joinpath(@__DIR__, "Makie", "CSG.jl"))
include(joinpath(@__DIR__, "Makie", "SolidStateDetector.jl"))
include(joinpath(@__DIR__, "Makie", "DriftPaths.jl"))
include(joinpath(@__DIR__, "Makie", "ChargeCloudModels.jl"))
include(joinpath(@__DIR__, "Makie", "Waveform.jl"))
include(joinpath(@__DIR__, "Makie", "Grids.jl"))

end # module
