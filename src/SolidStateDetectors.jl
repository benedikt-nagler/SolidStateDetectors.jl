# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

__precompile__(true)

module SolidStateDetectors

using LinearAlgebra
using Random
using Statistics

import Adapt
using ArraysOfArrays
using FillArrays
using Format
using Interpolations
using IntervalSets
using JSON
using KernelAbstractions
using LaTeXStrings
using LightXML
using ParallelProcessingTools
using ProgressMeter
using RadiationDetectorSignals
using RecipesBase
using Requires
using Rotations
using StaticArrays
using StatsBase
using Unitful
using UnitfulAtomic
using YAML

include("ka_compat.jl")

include("ConstructiveSolidGeometry/ConstructiveSolidGeometry.jl")
using .ConstructiveSolidGeometry
using .ConstructiveSolidGeometry:
            CartesianPoint, CartesianVector, CartesianZero, cartesian_zero, CylindricalPoint,
            AbstractCoordinatePoint, _convert_point,
            CartesianVector, AbstractCoordinateVector,
            Cartesian, Cylindrical, AbstractCoordinateSystem, CoordinateSystemType,
            Geometry, AbstractGeometry, AbstractSurfacePrimitive,
            parse_rotation_matrix, parse_translate_vector, parse_CSG_transformation,
            transform, CSG_dict, Transformations, combine_transformations, barycenter,
            LocalAffineFrame, cartesian_zero, global_frame, frame_transformation, 
            ConfigFileError, _parse_value, distance_squared, geom_round,
            LengthQuantity, AngleQuantity, get_scale

import .ConstructiveSolidGeometry: sample, to_internal_units, from_internal_units

export CartesianPoint, CartesianVector, CartesianZero, cartesian_zero, CylindricalPoint

import Clustering
import DataStructures
import Distributions
import SpecialFunctions
import GPUArrays
import IntervalSets
import OrderedCollections
import Tables
import TypedTables

import Base: size, sizeof, length, getindex, setindex!, axes, getproperty, broadcast,
             range, ndims, eachindex, enumerate, iterate, IndexStyle, eltype, in, convert,
             show, print, println, display, +, -, &, *

export SolidStateDetector
export SSD_examples

export Grid

export ElectricPotential, PointTypes, EffectiveChargeDensity, ImpurityScale, DielectricDistribution, WeightingPotential, ElectricField
export apply_initial_state!
export calculate_electric_potential!, calculate_weighting_potential!, calculate_electric_field!
export ElectricFieldChargeDriftModel, ADLChargeDriftModel, ADL2016ChargeDriftModel, IsotropicChargeDriftModel, InactiveLayerChargeDriftModel, PowerLawTemperatureModel
export ConstantImpurityDensity, LinearImpurityDensity, ThermalDiffusionLithiumDensity, PtypePNJunctionImpurityDensity
export LinBouleImpurityDensity, ParBouleImpurityDensity, LinExpBouleImpurityDensity, ParExpBouleImpurityDensity, SplineBouleImpurityDensity
export NoChargeTrappingModel, BoggsChargeTrappingModel, ConstantLifetimeChargeTrappingModel, CombinedChargeTrappingModel
export get_active_volume, is_depleted, estimate_depletion_voltage, set_point_type_depletion_handling!
export adjust_impurity_and_electric_potential_to_match_depletion!, adjust_bias_and_electric_potential!
export calculate_stored_energy, calculate_mutual_capacitance, calculate_capacitance_matrix
export simulate_waveforms
export run_geant4_simulation
export Simulation, simulate!
export Event, drift_charges!
export add_baseline_and_extend_tail
export NBodyChargeCloud

using Unitful: RealOrRealQuantity as RealQuantity
const SSDFloat = Union{Float16, Float32, Float64}

abstract type ChargeCarrier end
abstract type Electron <: ChargeCarrier end 
abstract type Hole <: ChargeCarrier end

include("examples.jl")

include("Units.jl")
include("Axes/DiscreteAxis.jl")
include("World/World.jl")
include("Grids/Grids.jl")

include("MaterialProperties/MaterialProperties.jl")
include("Config/Config.jl")
include("ChargeDensities/ChargeDensities.jl")
include("ImpurityDensities/ImpurityDensities.jl")
include("ChargeDriftModels/ChargeDriftModels.jl")

include("ScalarPotentials/ScalarPotential.jl")

include("ChargeTrapping/ChargeTrapping.jl")
include("SolidStateDetector/DetectorGeometries.jl")

include("PotentialCalculation/PotentialCalculation.jl")

include("ElectricField/ElectricField.jl")

include("ChargeCloudModels/ChargeCloudModels.jl")
include("ChargeDrift/ChargeDrift.jl")
include("SignalGeneration/SignalGeneration.jl")

include("ChargeStatistics/ChargeStatistics.jl")
include("ChargeClustering/ChargeClustering.jl")

include("Simulation/Simulation.jl")
include("Event/Event.jl")
include("MCEventsProcessing/MCEventsProcessing.jl")

include("IO/IO.jl")

include("PlotRecipes/PlotRecipes.jl")
export @P_str # protected strings to overwrite plot labels with units

# --- Makie plotting entry points ---------------------------------------------------------------
#
# The Makie recipes live in the `SolidStateDetectorsMakieExt` extension (weakdep `Makie`), but the
# plotting *functions* are declared here, so that there is exactly one function object per entry
# point: `SolidStateDetectors.electricfieldplot` **is** the recipe function, and the extension only
# attaches methods to it (`import SolidStateDetectors: electricfieldplot` before `@recipe`).
#
# A forwarding stub would be a different function object from the one Makie's `@recipe` generates,
# so `Makie.plot(obj)` / `plot!(ax, obj)` could never route to the recipe, and every method
# signature and docstring would be erased behind `args...; kwargs...`.
#
# Consequence: without a Makie backend loaded these functions have no methods at all. `__init__`
# registers a `MethodError` hint that says which package to load — see `_makie_backend_hint`.
export electricfieldplot, electricfieldplot!
function electricfieldplot end
function electricfieldplot! end

export scalarpotentialplot, scalarpotentialplot!
function scalarpotentialplot end
function scalarpotentialplot! end

export csgplot, csgplot!
function csgplot end
function csgplot! end

export detectorplot, detectorplot!
function detectorplot end
function detectorplot! end

export driftpathplot, driftpathplot!
function driftpathplot end
function driftpathplot! end

export chargecloudplot, chargecloudplot!
function chargecloudplot end
function chargecloudplot! end

export waveformplot, waveformplot!
function waveformplot end
function waveformplot! end

# The Makie counterpart of `plot_electron_and_hole_contribution`, under its own name rather than as
# a method on that function (D9): the Plots `@userplot` accepts `(args...; kwargs...)`, so a method
# specific enough to take `(::Event, ::Simulation, ::Int)` would capture every call the moment Makie
# is loaded, and the parity harness needs both backends usable in one session.
export ehcontributionplot, ehcontributionplot!
function ehcontributionplot end
function ehcontributionplot! end

export gridplot, gridplot!
function gridplot end
function gridplot! end

# The Makie counterpart of `plot_electric_fieldlines`, named for the same reason as
# `ehcontributionplot` (D9).
export fieldlinesplot, fieldlinesplot!
function fieldlinesplot end
function fieldlinesplot! end

const _makie_entry_points = (electricfieldplot, electricfieldplot!,
                             scalarpotentialplot, scalarpotentialplot!,
                             csgplot, csgplot!,
                             detectorplot, detectorplot!,
                             driftpathplot, driftpathplot!,
                             chargecloudplot, chargecloudplot!,
                             waveformplot, waveformplot!,
                             ehcontributionplot, ehcontributionplot!,
                             gridplot, gridplot!,
                             fieldlinesplot, fieldlinesplot!)

function _makie_backend_hint(io::IO, exc::MethodError, @nospecialize(argtypes), @nospecialize(kwargs))
    exc.f in _makie_entry_points || return nothing
    isnothing(Base.get_extension(@__MODULE__, :SolidStateDetectorsMakieExt)) || return nothing
    print(io, "\n\n`$(nameof(exc.f))` is provided by the SolidStateDetectors Makie extension. ",
        "Load a Makie backend first, e.g. `using CairoMakie` (or `using GLMakie`).")
    return nothing
end


export GridPoint
@deprecate GridPoint(grid::Grid{T, N, S}, idxs::Vararg{Int, N}) where {T, N, S} getpoint(grid, idxs...)
@deprecate GridPoint(grid::Grid{T, N, S}, idxs::NTuple{N, Int}) where {T, N, S} getpoint(grid, idxs)


function __init__()
    @require LegendHDF5IO ="c9265ca6-b027-5446-b1a4-febfa8dd10b0" begin
        include("../ext/SolidStateDetectorsLegendHDF5IOExt.jl")
    end
    Base.Experimental.register_error_hint(_makie_backend_hint, MethodError)
end

end # module
