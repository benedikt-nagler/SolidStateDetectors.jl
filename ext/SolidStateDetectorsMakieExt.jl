# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

module SolidStateDetectorsMakieExt

using Makie

using LinearAlgebra
using Unitful

using SolidStateDetectors
using SolidStateDetectors.ConstructiveSolidGeometry

# The plotting entry points are *declared* in the main module (see the "Makie plotting entry points"
# block in `src/SolidStateDetectors.jl`); `import` — not `using` — is what lets `@recipe` attach its
# methods to those functions instead of creating new ones in this module.
import SolidStateDetectors: electricfieldplot, electricfieldplot!

# Non-exported helpers / types reused from the main package.
import SolidStateDetectors: SSDFloat, ElectricField,
    get_crosssection_idx_and_value, get_extended_ticks, get_2π_potential,
    internal_efield_unit, internal_length_unit, internal_angle_unit, internal_voltage_unit
import SolidStateDetectors.ConstructiveSolidGeometry: Cartesian, Cylindrical
import StatsBase: midpoints

include(joinpath(@__DIR__, "Makie", "ElectricField.jl"))

end # module
