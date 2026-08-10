# Smoke tests for the Makie plotting extension (`SolidStateDetectorsMakieExt`).
#
# CairoMakie renders headless, so these run in CI without a display.
#
# NOTE: these are deliberately shallow — they cannot detect a wrong colormap, a missing colorbar, a
# transposed array or a lost title. As each recipe is ported they are replaced by assertions that
# compare title, axis labels, limits, `colorrange` and the plotted data against the output of the
# corresponding Plots.jl recipe.

using Test
using SolidStateDetectors
using CairoMakie

const ssd_makie_ext = Base.get_extension(SolidStateDetectors, :SolidStateDetectorsMakieExt)

@testset "Makie extension loaded" begin
    @test ssd_makie_ext !== nothing
    # The entry points are declared in the main module and given methods by the extension, so
    # they must be one and the same function (plan D4a) — a forwarding stub would not be.
    @test SolidStateDetectors.electricfieldplot === ssd_makie_ext.electricfieldplot
    @test !isempty(methods(SolidStateDetectors.electricfieldplot))
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

    # generic Makie entry points must route to the recipe (`plottype`/`convert_arguments`),
    # which is what keeps the `plot(sim.electric_field, ...)` call sites in the docs working
    gfap = Makie.plot(sim.electric_field, φ = 30)
    @test gfap isa Makie.FigureAxisPlot
    @test gfap.plot isa ssd_makie_ext.ElectricFieldPlot
    fig2 = Makie.Figure(); ax2 = Makie.Axis(fig2[1, 1])
    @test Makie.plot!(ax2, sim.electric_field, φ = 30) isa ssd_makie_ext.ElectricFieldPlot
end

@testset "ElectricField recipe (Cartesian)" begin
    simc = Simulation{Float32}(SSD_examples[:CGD])
    simulate!(simc, device_array_type = Array, refinement_limits = missing)
    @test electricfieldplot(simc.electric_field, x = 0) isa Makie.FigureAxisPlot
    @test electricfieldplot(simc.electric_field, y = 0) isa Makie.FigureAxisPlot
    @test electricfieldplot(simc.electric_field, z = 0) isa Makie.FigureAxisPlot
end
