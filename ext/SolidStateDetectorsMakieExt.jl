# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

module SolidStateDetectorsMakieExt

using Makie

using LinearAlgebra
using Unitful

using SolidStateDetectors
using SolidStateDetectors.ConstructiveSolidGeometry

# Non-exported helpers / types reused from the main package.
import SolidStateDetectors: SSDFloat, ElectricField,
    get_crosssection_idx_and_value, get_extended_ticks, get_2π_potential,
    internal_efield_unit, internal_length_unit, internal_angle_unit, internal_voltage_unit
import SolidStateDetectors.ConstructiveSolidGeometry: Cartesian, Cylindrical
import StatsBase: midpoints

include(joinpath(@__DIR__, "Makie", "ElectricField.jl"))

end # module
