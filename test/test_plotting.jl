# Smoke tests for the Makie plotting extension (`SolidStateDetectorsMakieExt`).
#
# CairoMakie renders headless, so these run in CI without a display. This file is only included
# (see `runtests.jl`) on Julia versions where the pinned Makie (0.21) can load — `using CairoMakie`
# must not be reached otherwise, as it would trigger a failing precompilation.

using Test
using SolidStateDetectors
using CairoMakie

@testset "Makie extension loaded" begin
    ext = Base.get_extension(SolidStateDetectors, :SolidStateDetectorsMakieExt)
    @test ext !== nothing
end

@testset "ElectricField recipe" begin
    sim = Simulation{Float32}(SSD_examples[:InvertedCoax])
    simulate!(sim, device_array_type = Array, refinement_limits = missing)

    # Cylindrical cross-sections: r/φ → Axis, z → PolarAxis
    fap = electricfieldplot(sim.electric_field, φ = 30)
    @test fap isa Makie.FigureAxisPlot
    @test fap.axis isa Makie.Axis
    @test electricfieldplot(sim.electric_field, φ = 30, full_det = true,
        contours_equal_potential = true) isa Makie.FigureAxisPlot
    @test electricfieldplot(sim.electric_field, r = 0.02) isa Makie.FigureAxisPlot
    fz = electricfieldplot(sim.electric_field, z = 0.04)
    @test fz isa Makie.FigureAxisPlot
    @test fz.axis isa Makie.PolarAxis

    # mutating form into an existing axis
    fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1])
    @test electricfieldplot!(ax, sim.electric_field, φ = 30) isa Makie.AbstractPlot
end

@testset "ElectricField recipe (Cartesian)" begin
    simc = Simulation{Float32}(SSD_examples[:CGD])
    simulate!(simc, device_array_type = Array, refinement_limits = missing)
    @test electricfieldplot(simc.electric_field, x = 0) isa Makie.FigureAxisPlot
    @test electricfieldplot(simc.electric_field, y = 0) isa Makie.FigureAxisPlot
    @test electricfieldplot(simc.electric_field, z = 0) isa Makie.FigureAxisPlot
end
