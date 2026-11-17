# Tests for the Makie plotting extension (`SolidStateDetectorsMakieExt`) and for the plot-parity
# harness it is gated by.
#
# CairoMakie renders headless, so these run in CI without a display.
#
# Two layers:
#
#   1. "parity harness" — the tier-0/tier-1 machinery itself: the committed golden specs load, the
#      spec round-trips through JSON, the block-mean digest actually catches a transposed array, and
#      the Makie introspector produces a comparable spec. These guard the gate; without them a
#      broken harness would silently pass every later step.
#   2. per-recipe tests. `ScalarPotentialPlot` is gated properly: every potentials case of the matrix
#      and `ElectricFieldPlot` are plotted, introspected and diffed against their golden specs by
#      `case_diffs` (see `parity/gate.jl`).

using Test
using SolidStateDetectors
using CairoMakie
using Unitful

include(joinpath(@__DIR__, "parity", "cases.jl"))
include(joinpath(@__DIR__, "parity", "makie_spec.jl"))
using .PlotParityCases
using .MakieSpec
using .MakieSpec.PlotParitySpec

const ssd_makie_ext = Base.get_extension(SolidStateDetectors, :SolidStateDetectorsMakieExt)

# The tier-1 comparison rules; needs `PlotParitySpec` and `Unitful` in scope.
include(joinpath(@__DIR__, "parity", "gate.jl"))

@testset "Makie extension loaded" begin
    @test ssd_makie_ext !== nothing
    # The entry points are declared in the main module and given methods by the extension, so
    # they must be one and the same function (plan D4a) — a forwarding stub would not be.
    for f in (:electricfieldplot, :scalarpotentialplot, :csgplot, :detectorplot,
              :driftpathplot, :chargecloudplot, :waveformplot, :ehcontributionplot, :gridplot,
              :fieldlinesplot)
        @test getfield(SolidStateDetectors, f) === getfield(ssd_makie_ext, f)
        @test !isempty(methods(getfield(SolidStateDetectors, f)))
        # …and the mutating half, which is what `plot!(ax, obj)` bottoms out at.
        bang = Symbol(f, :!)
        @test getfield(SolidStateDetectors, bang) === getfield(ssd_makie_ext, bang)
        @test !isempty(methods(getfield(SolidStateDetectors, bang)))
    end
end

@testset "parity harness" begin
    @testset "case matrix" begin
        cs = cases()
        @test length(cs) > 100
        @test length(unique(getfield.(cs, :id))) == length(cs)
        # Every case id must be usable as a filename.
        @test all(id -> occursin(r"^[A-Za-z0-9_]+$", id), getfield.(cs, :id))
        # `r` is a Plots alias for `:rotation`, so those cases carry a note and must never be used
        # as a parity reference (see PlotParityCases.NOTE_R_ALIAS).
        @test all(c -> !haskey(c.kwargs, :r) || !isempty(c.note), cs)
        @test !isempty(cases(groups = ["potentials"]))
        @test all(c -> !c.slow, cases(include_slow = false))
    end

    @testset "golden specs" begin
        cs = cases()
        missing_specs = [c.id for c in cs if !PlotParitySpec.has_spec(c.id)]
        @test isempty(missing_specs)

        spec = read_spec("potentials_cyl_electricpotential_phi")
        @test spec.backend == "plots"
        @test isempty(spec.error)
        @test length(spec.axes) == 1
        ax = spec.axes[1]
        # The properties a Makie recipe cannot set, and which therefore define the work of the
        # figure-level wrappers (plan §1.1).
        @test ax.title == "Electric Potential @ φ = 0.0°"
        @test (ax.xguide, ax.yguide) == ("r", "z")
        @test (ax.xunit, ax.yunit, ax.zunit) == ("m", "m", "V")
        @test ax.clims ≈ [0.0, 3500.0015]   # rounded to 8 significant digits on write
        @test ax.aspect_ratio == "1"
        @test length(spec.series) == 1
        @test spec.series[1].seriestype == "heatmap"

        # The cylindrical z cross-section is the one that must become a PolarAxis.
        @test read_spec("potentials_cyl_electricpotential_z").axes[1].projection == "polar"
    end

    @testset "spec round trip" begin
        dir = mktempdir()
        original = read_spec("efield_cyl_phi")
        write_spec(original; dir)
        roundtripped = read_spec("efield_cyl_phi"; dir)
        @test roundtripped.id == original.id
        @test roundtripped.axes[1].title == original.axes[1].title
        @test roundtripped.axes[1].clims == original.axes[1].clims
        @test roundtripped.series[1].z.size == original.series[1].z.size
        @test roundtripped.series[1].z.digest ≈ original.series[1].z.digest
    end

    @testset "digest catches a transposed array" begin
        # A non-symmetric matrix: transposing it leaves size (square), extrema and even the sorted
        # values untouched, so only the block-mean digest can tell the two apart. This is the check
        # that makes plan §1.7 testable rather than a comment.
        a = Float64[i + 10j for i in 1:24, j in 1:24]
        upright, flipped = data_spec(a), data_spec(permutedims(a))
        @test upright.size == flipped.size
        @test upright.extrema == flipped.extrema
        @test !isempty(compare_data(upright, flipped))
        @test isempty(compare_data(upright, flipped; transposed = true))
        @test isempty(compare_data(upright, upright))

        # Non-square arrays are caught on shape alone, and `transposed` reconciles them.
        b = Float64[i + 10j for i in 1:20, j in 1:24]
        @test !isempty(compare_data(data_spec(b), data_spec(permutedims(b))))
        @test isempty(compare_data(data_spec(b), data_spec(permutedims(b)); transposed = true))

        # Data that cannot be reduced to numbers is recorded as absent rather than recursed into:
        # `Float64.(::Vector{Point{3}})` is still a vector of points, which used to loop forever.
        @test data_spec([Makie.Point3f(1, 2, 3), Makie.Point3f(4, 5, 6)]) === nothing
    end

    @testset "point-based primitives are split into x/y/z" begin
        # A Makie line carries one `Vector{Point}` where Plots carried separate x/y/z vectors, so
        # the introspector splits it — without that, no line-based case (drift paths, waveforms, CSG)
        # could ever be compared against its golden spec.
        fig = Makie.Figure()
        ax = Makie.Axis(fig[1, 1])
        Makie.lines!(ax, [1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
        spec = makie_case_spec(fig; id = "probe", group = "probe", call = "probe")
        s = only(spec.series)
        @test s.x.extrema == [1.0, 3.0]
        @test s.y.extrema == [4.0, 6.0]
        @test s.z === nothing
    end

    @testset "title conversion" begin
        # The gate converts the golden title's cut position into display units before comparing,
        # because that is the one deliberate difference in the string (the axes are in mm, so the
        # title is too). Everything else about it still has to match character for character.
        @test _expected_title("Electric Potential @ z = 0.0375 m"; ext = ssd_makie_ext) ==
              "Electric Potential @ z = 37.5 mm"
        @test _expected_title("Point Type Map @ φ = 0.0°"; ext = ssd_makie_ext) ==
              "Point Type Map @ φ = 0.0°"                      # degrees are already display units
        @test _expected_title("Dielectric Distribution"; ext = ssd_makie_ext) ==
              "Dielectric Distribution"                        # no position to convert
    end

    @testset "colormap comparison" begin
        viridis = "grad(#440154,#21908c,#fde725)"
        @test isempty(compare_colormap(viridis, viridis))
        # A stop or two of interpolation difference between the backends is tolerated…
        @test isempty(compare_colormap(viridis, "grad(#440154,#219a8e,#fde725)"))
        # …a different colormap is not.
        @test !isempty(compare_colormap(viridis, "grad(#000004,#bb3755,#fcffa4)"))
    end

    @testset "extent check catches a rescaled geometry" begin
        # The only data comparison a geometry case gets: no heatmap to diff, and a mesh is held so
        # differently by the two backends that neither the series count nor any single array is
        # comparable (see `_extent_diffs`). It has to be non-vacuous, and it has to bridge the units.
        box = SolidStateDetectors.ConstructiveSolidGeometry.Box{Float32}(hX = 1, hY = 2, hZ = 1)
        golden = read_spec("csg_box")
        actual = makie_case_spec(Makie.plot(box); id = "csg_box", group = "csg", call = "probe")
        # Golden metres against Makie millimetres: the check is only meaningful because it converts.
        @test _data_bbox(golden, 1) == [-1.0, 1.0]
        @test _data_bbox(actual, 1) == [-1000.0, 1000.0]
        @test isempty(_extent_diffs(golden, actual; ext = ssd_makie_ext))

        flat = SolidStateDetectors.ConstructiveSolidGeometry.Box{Float32}(hX = 1, hY = 2, hZ = 0.5)
        shrunk = makie_case_spec(Makie.plot(flat); id = "csg_box", group = "csg", call = "probe")
        @test !isempty(_extent_diffs(golden, shrunk; ext = ssd_makie_ext))

        # …and `case_diffs` must actually *reach* it. It once did not: the `_extent_diffs` branch sat
        # behind an earlier `if gh === nothing || ah === nothing`, so it was unreachable and every
        # geometry case was gated on its title and labels alone. A unit test of the function cannot
        # see that; only driving it through the gate can.
        reached, _ = case_diffs(golden, shrunk; ext = ssd_makie_ext)
        @test any(d -> occursin("extent", d), reached)
        @test isempty(first(case_diffs(golden, actual; ext = ssd_makie_ext)))
    end

    @testset "colormaps are recorded only where they are used" begin
        # Every Makie plot type carries a `colormap` attribute, viridis by default. Reading it off a
        # flat-coloured `mesh!` gave every geometry case a colormap the picture never uses, and the
        # gate then failed all 31 of them against a golden that correctly recorded none.
        box = SolidStateDetectors.ConstructiveSolidGeometry.Box{Float32}(hX = 1, hY = 1, hZ = 1)
        @test makie_case_spec(Makie.plot(box); id = "probe", group = "csg",
                              call = "probe").axes[1].colormap == ""
    end
end

@testset "shared infrastructure" begin
    ext = ssd_makie_ext

    @testset "display units and labels (D2)" begin
        # Lengths read better in mm than in the internal metres; angles in degrees; waveform times
        # in ns. Everything else stays in its internal unit.
        @test ext.display_unit(u"m") == u"mm"
        @test ext.display_unit(u"mm") == u"mm"
        @test ext.display_unit(u"rad") == u"°"
        @test ext.display_unit(u"°") == u"°"
        @test ext.display_unit(u"s") == u"ns"
        @test ext.display_unit(u"V") == u"V"
        @test ext.display_unit(u"V/m") == u"V/m"
        @test ext.display_unit(Unitful.NoUnits) == Unitful.NoUnits
        @test ext.display_unit(nothing) == Unitful.NoUnits

        # Limits and other axis-space numbers must be converted, because an axis carrying a
        # UnitfulConversion works in display space.
        @test ext.display_value(0.04u"m") ≈ 40.0
        @test ext.display_value(Float32(π)u"rad") ≈ 180.0 rtol = 1e-6
        @test ext.display_value(3.0) ≈ 3.0

        # `unitformat = :slash` reproduced exactly; `:parens` is the form Plots uses for colorbars.
        @test ext.unit_label("r", u"mm") == "r / mm"
        @test ext.unit_label("z", u"m") == "z / m"
        @test ext.unit_label("Weighting Potential", Unitful.NoUnits) == "Weighting Potential"
        @test ext.unit_label("Electric Potential", u"V"; format = :parens) == "Electric Potential (V)"
        @test_throws ArgumentError ext.unit_label("r", u"mm"; format = :nonsense)

        @test ext.unit_conversion(Unitful.NoUnits) === nothing   # Makie's own default
        @test ext.unit_conversion(u"m") isa Makie.UnitfulConversion

        # Heatmap data is rescaled by hand rather than by an axis conversion, because Makie 0.24's
        # `heatmap` cannot take Unitful data at all — the two calls below are what would break.
        @test ext.display_factor(u"m") ≈ 1000
        @test ext.display_factor(u"rad") ≈ 180 / π
        @test ext.display_factor(u"V") == 1
        @test ext.display_values(Float32[0, 0.04], u"m") ≈ [0.0, 40.0]
        @test eltype(ext.display_values(Float32[0, 0.04], u"m")) === Float32
        @test ext.display_values([0, 1], Unitful.NoUnits) == [0, 1]
        # (The failure mode itself is not exercised here: `heatmap!(ax, xs*u"m", …)` overflows the
        # stack rather than throwing, which is not something to provoke inside a test process.)
    end

    @testset "axis construction" begin
        # The properties a recipe cannot set and a wrapper must: title, labels, aspect, limits —
        # and the limits must come out in display space (mm), not in the data's metres.
        fig = Makie.Figure()
        ax = ext.ssd_axis(fig[1, 1]; title = "Electric Potential @ φ = 0.0°",
                          xname = "r", yname = "z",
                          xunit = u"m", yunit = u"m",
                          xlims = (0.0u"m", 0.04u"m"), ylims = (-0.01u"m", 0.09u"m"),
                          aspect_ratio_one = true)
        @test ax isa Makie.Axis
        @test ax.title[] == "Electric Potential @ φ = 0.0°"
        @test (ax.xlabel[], ax.ylabel[]) == ("r / mm", "z / mm")
        @test ax.aspect[] isa Makie.DataAspect
        Makie.update_state_before_display!(fig)
        lims = ax.finallimits[]
        @test lims.origin[1] ≈ 0.0 atol = 1e-6
        @test lims.origin[1] + lims.widths[1] ≈ 40.0 atol = 1e-6
        @test lims.origin[2] ≈ -10.0 atol = 1e-6
        @test lims.origin[2] + lims.widths[2] ≈ 90.0 atol = 1e-6

        # Bare numbers are interpreted as being in the data's unit, so they land in the same place.
        fig2 = Makie.Figure()
        ax2 = ext.ssd_axis(fig2[1, 1]; xname = "r", xunit = u"m", xlims = (0.0, 0.04))
        Makie.update_state_before_display!(fig2)
        @test ax2.finallimits[].widths[1] ≈ 40.0 atol = 1e-6

        # Unit conversions are *off* by default: the recipes hand the axis numbers that are already
        # in display units, and an axis carrying a conversion rejects those (see `display_values`).
        @test ax.dim1_conversion[] === nothing || ax.dim1_conversion[] isa Makie.NoDimConversion
        axc = ext.ssd_axis(Makie.Figure()[1, 1]; xname = "r", xunit = u"m", dim_conversions = true)
        @test axc.dim1_conversion[] isa Makie.UnitfulConversion

        # A PolarAxis has no labels, no aspect and no Cartesian limits; asking for them must not
        # throw, because the cylindrical z cross-section goes through exactly this path.
        figp = Makie.Figure()
        axp = ext.ssd_axis(figp[1, 1]; axistype = Makie.PolarAxis, title = "Point Type Map @ z = 0",
                           xname = "φ", yname = "r", xunit = u"rad", yunit = u"m",
                           aspect_ratio_one = true)
        @test axp isa Makie.PolarAxis
        @test axp.title[] == "Point Type Map @ z = 0"

        # An Axis3 carries the third label, equal data scaling as the symbol `:data`, and
        # `viewmode = :fit`. The last one is not cosmetic: under Makie's default `:fitzoom` an
        # elongated geometry pushes the z label clean off the figure, with the attribute still set
        # correctly — invisible to the gate, and caught by the tier-2 gallery.
        fig3 = Makie.Figure()
        ax3 = ext.ssd_axis(fig3[1, 1]; axistype = Makie.Axis3, xname = "x", yname = "y",
                           zname = "z", xunit = u"m", yunit = u"m", zunit = u"m",
                           aspect_ratio_one = true)
        @test ax3 isa Makie.Axis3
        @test (ax3.xlabel[], ax3.ylabel[], ax3.zlabel[]) == ("x / mm", "y / mm", "z / mm")
        @test ax3.aspect[] === :data
        @test ax3.viewmode[] === :fit

        # An Axis3 takes six limits, not four, so the third dimension has to reach it.
        ax3l = ext.ssd_axis(Makie.Figure()[1, 1]; axistype = Makie.Axis3,
                            xunit = u"m", yunit = u"m", zunit = u"m",
                            xlims = (0.0, 0.04), ylims = (0.0, 0.04), zlims = (-0.01, 0.09))
        @test ax3l.limits[] == (0.0, 40.0, 0.0, 40.0, -10.0, 90.0)
    end

    @testset "colorbar" begin
        # A Colorbar built before display would read `Automatic()` off the plot, so the explicit
        # colorrange is what makes it correct (plan §1.4).
        fig = Makie.Figure()
        ax = Makie.Axis(fig[1, 1])
        hm = Makie.heatmap!(ax, rand(Float32, 8, 6))
        cb = ext.add_colorbar!(fig, hm; label = "Electric Potential / V", colorrange = (0.0, 3500.0))
        @test cb isa Makie.Colorbar
        @test cb.label[] == "Electric Potential / V"
        @test cb.colorrange[] == (0.0, 3500.0)
        @test any(b -> b isa Makie.Colorbar, fig.content)

        # The wrappers resolve `Makie.automatic` themselves (the bar needs a range before display),
        # so their resolution has to be Makie's own — including on constant data, where an explicit
        # equal-ended colorrange throws while `automatic` renders.
        @test ext.color_limits([1.0 2.0; 3.0 4.0]) == (1.0, 4.0)
        @test ext.color_limits((0.0, 3500.0)) == (0.0, 3500.0)
        @test ext.color_limits(zeros(Float32, 3, 4)) == (-0.5, 0.5)   # Makie's distinct_extrema_nan
        constant = Makie.Figure()
        constant_ax = Makie.Axis(constant[1, 1])
        flat = Makie.heatmap!(constant_ax, zeros(Float32, 4, 3);
                              colorrange = ext.color_limits(zeros(Float32, 4, 3)))
        @test ext.add_colorbar!(constant, flat; colorrange = ext.color_limits(zeros(Float32, 4, 3))) isa Makie.Colorbar
        @test Makie.colorbuffer(constant) isa AbstractMatrix   # renders rather than throwing

        # A composite plot carrying a heatmap *and* a contour overlay has two colormaps, and
        # `Colorbar` refuses to guess; the colour-mapped child is what it must be built from.
        @test ext.colorbar_source(hm) === hm
        comp = scalarpotentialplot(fixtures().cyl.electric_potential;
                                   φ = 0, contours_equal_potential = true)
        @test Makie.plotkey(ext.colorbar_source(comp.plot)) === :heatmap
        @test count(b -> b isa Makie.Colorbar, comp.figure.content) == 1
    end

    @testset "palette (D3)" begin
        # Pinned to Plots' default palette so contacts keep the colours the golden specs recorded.
        @test ext.contact_color(1) == Makie.to_color("#009afa")
        @test ext.contact_color(2) == Makie.to_color("#e36f47")
        @test ext.contact_color(3) == Makie.to_color("#3ea44e")
        # …and demonstrably *not* Makie's default cycle, which is what a naive port would have used.
        @test ext.contact_color(1) != Makie.wong_colors()[1]
        # Wraps around rather than erroring on a detector with many contacts.
        @test ext.contact_color(length(ext.SSD_PALETTE) + 1) == ext.contact_color(1)
        @test ext.ssd_theme() isa Makie.Theme

        golden = read_spec("detector_cyl_contact1")
        @test golden.series[1].color == "#009afa"
    end

    @testset "legend elements" begin
        sim = Simulation{Float32}(SSD_examples[:InvertedCoax])
        det = sim.detector

        elements, labels = ext.legend_entries(det)
        # Same entries, in the same order, with the same text as the Plots legend.
        @test labels == read_spec("detector_cyl_csg").axes[1].legend_entries
        @test length(elements) == length(labels)
        @test all(e -> e isa Makie.PolyElement, elements)

        _, with_semi = ext.legend_entries(det; show_semiconductor = true)
        @test "Semiconductor" in with_semi
        _, no_passives = ext.legend_entries(det; show_passives = false)
        @test !any(l -> occursin("Passiv", l), no_passives)

        de, dl = ext.drift_legend_entries()
        @test dl == read_spec("driftpaths_single").axes[1].legend_entries == ["e-", "h+"]
        @test all(e -> e isa Makie.LineElement, de)

        fig = Makie.Figure()
        Makie.Axis(fig[1, 1])
        @test ext.add_legend!(fig, elements, labels) isa Makie.Legend
        @test ext.add_legend!(fig, Makie.LegendElement[], String[]) === nothing
    end
end

@testset "Potentials" begin
    ext = ssd_makie_ext
    potential_cases = cases(groups = ["potentials"])

    @testset "tier-1 parity gate" begin
        # Every potentials case of the matrix, plotted, introspected and diffed against the golden
        # spec extracted from Plots. This is the gate of step 3: it fails on a lost title, a wrong
        # colormap, a missing colorbar, a shifted limit or a transposed array — none of which a
        # smoke test can see.
        gated = 0
        for c in potential_cases
            golden = read_spec(c.id)
            # `r =` cross sections have no valid Plots reference (`PlotParityCases.NOTE_R_ALIAS`);
            # they are asserted directly below instead.
            isempty(golden.note) && isempty(golden.error) || continue
            args = c.args()
            fig = scalarpotentialplot(args...; c.kwargs...)
            spec = makie_case_spec(fig; id = c.id, group = c.group, call = call_string(c, args))
            diffs, _ = case_diffs(golden, spec; ext)
            @test (c.id, diffs) == (c.id, String[])
            gated += 1
        end
        @test gated == 41   # 49 cases minus the 8 that `r` aliasing makes unusable as a reference
    end

    @testset "r cross sections" begin
        # Unreachable through Plots (`r` is its alias for `:rotation`), so there is no golden spec —
        # these are asserted against the data instead. Makie has no such alias, and this is the first
        # backend on which an `r =` cross section works at all.
        sim = fixtures().cyl
        sp = sim.electric_potential
        fap = scalarpotentialplot(sp; r = 20u"mm")
        spec = makie_case_spec(fap; id = "potentials_cyl_electricpotential_r_units",
                               group = "potentials", call = "probe")
        ax = spec.axes[1]
        # 20 mm is not a grid tick; like every other cross section it snaps to the closest one
        # (19.355 mm, which the title rounds to three significant digits).
        @test ax.title == "Electric Potential @ r = 19.4 mm"
        @test (ax.xguide, ax.yguide) == ("φ / °", "z / mm")
        @test ax.projection == "none"
        @test ax.colorbar_title == "Electric Potential (V)"

        # The drawn slice is the one at the requested radius, in the (φ, z) orientation Makie wants.
        expanded = ext._full_φ(sp)
        _, idx, value = SolidStateDetectors.get_crosssection_idx_and_value(
            expanded.grid, 0.02, missing, missing)
        @test value ≈ 0.0194 atol = 1e-4
        heat = only(filter(s -> s.seriestype == "heatmap", spec.series))
        @test isempty(compare_data(data_spec(expanded.data[idx, :, :]), heat.z, label = "z"))

        # An `r =` cross section is the cylindrical *shell* at that radius, unrolled into (φ, z) —
        # not a cut through the detector. Its axis box therefore carries the shell's proportions,
        # 2πr wide by Δz tall (the fixture spans z = -10…90 mm), so the picture is not stretched
        # across the full figure width.
        shell = only(filter(b -> b isa Makie.Axis, fap.figure.content))
        @test shell.aspect[] isa Makie.AxisAspect
        @test shell.aspect[].aspect ≈ 2π * ext.display_values(value, u"m") / 100 rtol = 1e-3

        # The title is in display units whatever the caller passed, so it agrees with the axis
        # labels; `0.02` (internal metres) and `20u"mm"` give the same string.
        @test makie_case_spec(scalarpotentialplot(sp; r = 0.02); id = "probe", group = "potentials",
                              call = "probe").axes[1].title == "Electric Potential @ r = 19.4 mm"
    end

    @testset "entry points" begin
        sp = fixtures().cyl.electric_potential
        # `plot(sp)` must keep producing what Plots produced: title, labels and colorbar, not a bare
        # heatmap — that is what the doc call sites rely on.
        fap = Makie.plot(sp, φ = 0)
        @test fap isa Makie.FigureAxisPlot
        @test fap.plot isa ssd_makie_ext.ScalarPotentialPlot
        @test any(b -> b isa Makie.Colorbar, fap.figure.content)
        @test fap.axis.title[] == "Electric Potential @ φ = 0.0°"

        # …while the mutating form belongs to the caller's axis and sets nothing on it.
        fig = Makie.Figure()
        ax = Makie.Axis(fig[1, 1])
        plt = scalarpotentialplot!(ax, sp; φ = 0)
        @test plt isa ssd_makie_ext.ScalarPotentialPlot
        @test ax.title[] == ""
        @test Makie.plot!(ax, sp, φ = 0) isa ssd_makie_ext.ScalarPotentialPlot

        # The cylindrical z cross section is the one that has to become a PolarAxis.
        zfap = scalarpotentialplot(sp; z = 0.04)
        @test zfap.axis isa Makie.PolarAxis
        @test Makie.preferred_axis_type(Makie.plot(sp, z = 0.04).plot) === Makie.PolarAxis
        # …and it must be drawn as polar cells: a heatmap's axis-aligned cells collapse into a
        # sliver under the polar transform, which the tier-1 gate cannot see (it compares the data,
        # not the geometry) and the tier-2 gallery caught.
        zspec = makie_case_spec(zfap; id = "probe", group = "potentials", call = "probe")
        @test _field_series(zspec).seriestype == "polarcells"

        # The exact data check, against the fixture: the drawn cells are the slice, in (φ, r) order,
        # one row per φ tick, with nothing interpolated between them.
        cells = ssd_makie_ext.colorbar_source(zfap.plot)
        drawn = MakieSpec.polar_cell_matrix(cells)
        expanded = ext._full_φ(sp)
        _, zidx, _ = SolidStateDetectors.get_crosssection_idx_and_value(
            expanded.grid, missing, missing, 0.04)
        @test drawn == permutedims(expanded.data[:, :, zidx])

        # Attributes still win over the field's own defaults.
        custom = scalarpotentialplot(sp; φ = 0, colormap = :turbo, colorrange = (0, 1000))
        hm = ssd_makie_ext.colorbar_source(custom.plot)
        @test hm.colormap[] === :turbo
        @test collect(hm.colorrange[]) ≈ [0.0, 1000.0]
        # …and the colorbar has to describe the same scale as the plot it sits next to.
        cb = only(filter(b -> b isa Makie.Colorbar, custom.figure.content))
        @test collect(cb.colorrange[]) ≈ [0.0, 1000.0]
    end
end

@testset "Electric field magnitude" begin
    ext = ssd_makie_ext
    field_cases = cases(groups = ["efield"])

    @testset "tier-1 parity gate" begin
        gated = 0
        for c in field_cases
            golden = read_spec(c.id)
            # Plots intercepts `r` as its rotation alias, so the two radial cases have no valid
            # reference. The equivalent scalar-potential cases are checked directly above.
            isempty(golden.note) && isempty(golden.error) || continue
            args = c.args()
            fig = electricfieldplot(args...; c.kwargs...)
            spec = makie_case_spec(fig; id = c.id, group = c.group, call = call_string(c, args))
            diffs, _ = case_diffs(golden, spec; ext)
            @test (c.id, diffs) == (c.id, String[])
            gated += 1
        end
        @test gated == 10   # 12 cases minus the 2 radial cases affected by Plots' `r` alias
    end

    @testset "entry points" begin
        ef = fixtures().cyl.electric_field
        fap = Makie.plot(ef, φ = 0)
        @test fap isa Makie.FigureAxisPlot
        @test fap.plot isa ext.ElectricFieldPlot
        @test fap.axis.title[] == "Electric Field (Magn.) @ φ = 0.0°"
        @test any(b -> b isa Makie.Colorbar, fap.figure.content)

        fig = Makie.Figure()
        ax = Makie.Axis(fig[1, 1])
        plt = electricfieldplot!(ax, ef; φ = 0)
        @test plt isa ext.ElectricFieldPlot
        @test ax.title[] == ""
        @test Makie.plot!(ax, ef, φ = 0) isa ext.ElectricFieldPlot

        zfap = electricfieldplot(ef; z = 0.04)
        @test zfap.axis isa Makie.PolarAxis
        @test Makie.preferred_axis_type(Makie.plot(ef, z = 0.04).plot) === Makie.PolarAxis
        zspec = makie_case_spec(zfap; id = "probe", group = "efield", call = "probe")
        @test _field_series(zspec).seriestype == "polarcells"

        custom = electricfieldplot(ef; φ = 0, colormap = :turbo, colorrange = (0, 1e6))
        hm = ext.colorbar_source(custom.plot)
        @test hm.colormap[] === :turbo
        @test collect(hm.colorrange[]) ≈ [0.0, 1e6]
        cb = only(filter(b -> b isa Makie.Colorbar, custom.figure.content))
        @test collect(cb.colorrange[]) ≈ [0.0, 1e6]
    end

    @testset "colour scale describes the drawn slice" begin
        # The Plots recipe sets no `clims` for a field magnitude, so GR scaled every panel to the
        # slice it drew. A Makie `Colorbar` must be handed a range explicitly (plan §1.4), and the
        # tier-1 gate cannot check the one the heatmap uses — `plot.colorrange[]` stays
        # `Automatic()` even after layout, so the spec records no colorrange at all. Hence this
        # assertion: the bar, the heatmap and the plotted data must all be on one scale.
        car = fixtures().car.electric_field
        drawn_maxima = Float64[]
        for kw in ((; x = 0), (; y = 0), (; z = 0))
            fap = electricfieldplot(car; kw...)
            hm = ext.colorbar_source(fap.plot)
            cb = only(filter(b -> b isa Makie.Colorbar, fap.figure.content))
            drawn = collect(extrema(hm[3][]))          # the matrix actually handed to the heatmap
            @test collect(hm.colorrange[]) ≈ drawn
            @test collect(cb.colorrange[]) ≈ drawn
            push!(drawn_maxima, drawn[2])
        end
        # …and the distinction is not academic: on this fixture the widest cross section still
        # covers well under the whole field's range, so a colorbar built from `extrema(norm.(ef))`
        # would advertise a scale roughly three times the one the picture uses.
        whole_max = maximum(sqrt(sum(abs2, v)) for v in car.data)
        @test maximum(drawn_maxima) < whole_max

        # The cylindrical branches go through the same helper, including the mirrored `full_det`
        # slice, which is the one place where the drawn array is not a plain view of the field.
        cyl = fixtures().cyl.electric_field
        for kw in ((; φ = 0), (; φ = 0, full_det = true), (; z = 0.04))
            fap = electricfieldplot(cyl; kw...)
            plt = ext.colorbar_source(fap.plot)
            cb = only(filter(b -> b isa Makie.Colorbar, fap.figure.content))
            drawn = collect(extrema(Makie.plotkey(plt) === :mesh ? plt.color[] : plt[3][]))
            @test collect(plt.colorrange[]) ≈ drawn
            @test collect(cb.colorrange[]) ≈ drawn
        end
    end
end

@testset "CSG" begin
    ext = ssd_makie_ext
    CSG = SolidStateDetectors.ConstructiveSolidGeometry

    box  = CSG.Box{Float32}(hX = 1, hY = 2, hZ = 1)
    cone = CSG.Cone{Float32}(r = 1, hZ = 1)
    leaves(fap) = MakieSpec.leaf_plots(fap isa Makie.FigureAxisPlot ? fap.plot : fap)
    keys_of(fap) = [Makie.plotkey(l) for l in leaves(fap)]
    only_of(fap, key) = only(filter(l -> Makie.plotkey(l) === key, leaves(fap)))

    @testset "tier-1 parity gate" begin
        gated = 0
        for c in cases(groups = ["csg"])
            golden = read_spec(c.id)
            # `plot(::CartesianVector)` errors under Plots, so that case has no reference
            # (`NOTE_VECTOR_BROKEN`); it is asserted directly below.
            isempty(golden.note) && isempty(golden.error) || continue
            args = c.args()
            fig = csgplot(args...; c.kwargs...)
            spec = makie_case_spec(fig; id = c.id, group = c.group, call = call_string(c, args))
            diffs, _ = case_diffs(golden, spec; ext)
            @test (c.id, diffs) == (c.id, String[])
            gated += 1
        end
        @test gated == 30   # 31 cases minus `csg_vector`, which Plots cannot draw at all
    end

    @testset "representations" begin
        # What each `seriestype` draws. The gate compares the extent of the geometry but not how it
        # is built, so a `:wireframe` that quietly also filled its faces would pass it.
        @test keys_of(csgplot(box)) == [:mesh, :linesegments]
        @test keys_of(csgplot(box; seriestype = :wireframe)) == [:linesegments]
        @test keys_of(csgplot(box; seriestype = :mesh3d)) == [:mesh, :linesegments]
        @test keys_of(csgplot(box; seriestype = :samplesurface)) == [:scatter]
        @test keys_of(csgplot(box; seriestype = :slice, z = 0)) == [:scatter]
        @test_throws ArgumentError csgplot(box; seriestype = :nonsense)

        # `:csg` puts a black wireframe over unlit faces, `:mesh3d` the thinner cell grid instead:
        # the two differ only in the line layer. The grid is black as well. The Plots recipe asked
        # for `linecolor --> :white`, but GR drew it dark, and a white grid on translucent faces
        # over a white background reads as a halo around the mesh rather than as a grid.
        @test only_of(csgplot(box), :linesegments).color[] == Makie.to_color(:black)
        mesh_grid = only_of(csgplot(box; seriestype = :mesh3d), :linesegments)
        @test mesh_grid.color[] == Makie.to_color(:black)
        @test mesh_grid.linewidth[] == ext.SSD_LINEWIDTHS.meshgrid

        # A box is six quadrangles; every face is fanned into triangles, since Makie meshes are
        # triangular and SSD's faces are quads, triangles or (on a prism) N-gons.
        m = only_of(csgplot(box), :mesh)[1][]
        @test length(Makie.GeometryBasics.coordinates(m)) == 24     # 4 vertices × 6 faces
        @test length(Makie.GeometryBasics.faces(m)) == 12           # 2 triangles × 6 faces
        @test all(f -> length(f) == 3, Makie.GeometryBasics.faces(m))
    end

    @testset "a subtracted primitive is drawn as a shell" begin
        # The one thing neither the gate's extent check nor a label can see: `box - cone` has the
        # same bounding box as `box`, and differs only in that the cone is drawn as a translucent
        # white shell with a thinner outline (`fillcolor := :white` in the Plots recipe).
        fap = csgplot(box - cone)
        meshes = filter(l -> Makie.plotkey(l) === :mesh, leaves(fap))
        @test length(meshes) == 2
        @test meshes[1].color[] == ext.contact_color(1)
        @test meshes[1].alpha[] == ext.SSD_ALPHAS.mesh
        @test meshes[2].color[] == Makie.to_color(:white)
        @test meshes[2].alpha[] == ext.SSD_ALPHAS.contact

        wires = filter(l -> Makie.plotkey(l) === :linesegments, leaves(fap))
        @test length(wires) == 2
        @test wires[2].linewidth[] ≈ ext.SSD_LINEWIDTHS.subtracted   # stored as Float32

        # A plain union has no subtracted part, so it draws one mesh for both primitives.
        @test length(filter(l -> Makie.plotkey(l) === :mesh, leaves(csgplot(box + cone)))) == 1
    end

    @testset "slices" begin
        # A `φ` cut is taken by rotating the geometry and cutting at `y = 0`, then keeping the half
        # plane `r >= 0` — the negative half is the φ + 180° cross section, and `full_det` is what
        # asks for both.
        half = only_of(csgplot(cone; seriestype = :slice, φ = 0), :scatter)[1][]
        full = only_of(csgplot(cone; seriestype = :slice, φ = 0, full_det = true), :scatter)[1][]
        @test minimum(p -> p[1], half) >= 0
        @test minimum(p -> p[1], full) < 0
        @test length(full) > length(half)

        # The axes are named after the cut, not after the Cartesian coordinates the points are
        # sampled in: a `φ` cut is drawn in (x, z) and labelled (r, z).
        @test (csgplot(cone; seriestype = :slice, φ = 0).axis.xlabel[],
               csgplot(cone; seriestype = :slice, φ = 0).axis.ylabel[]) == ("r / mm", "z / mm")
        @test csgplot(box; seriestype = :slice, z = 0).axis.xlabel[] == "x / mm"
        @test csgplot(box; seriestype = :slice, z = 0).axis.ylabel[] == "y / mm"
        @test csgplot(box; seriestype = :slice, x = 0).axis.ylabel[] == "z / mm"
    end

    @testset "entry points" begin
        # `plot(geom)` must keep producing the labelled figure with its legend, not a bare mesh.
        fap = Makie.plot(box)
        @test fap isa Makie.FigureAxisPlot
        @test fap.plot isa ext.CSGPlot
        @test fap.axis isa Makie.Axis3
        @test (fap.axis.xlabel[], fap.axis.ylabel[], fap.axis.zlabel[]) ==
              ("x / mm", "y / mm", "z / mm")
        # `aspect = :data`, so a 1 × 2 × 1 box is drawn 1 × 2 × 1 rather than in GR's fixed cube.
        @test fap.axis.aspect[] === :data

        # …while the mutating form belongs to the caller's axis and sets nothing on it.
        fig = Makie.Figure()
        ax = Makie.Axis3(fig[1, 1])
        @test csgplot!(ax, box) isa ext.CSGPlot
        @test ax.title[] == ""
        @test Makie.plot!(ax, box) isa ext.CSGPlot

        # A slice is the only representation that is not geometry in space.
        @test Makie.preferred_axis_type(Makie.plot(box, seriestype = :slice, z = 0).plot) ===
              Makie.Axis
        @test Makie.preferred_axis_type(Makie.plot(box).plot) === Makie.Axis3

        # Every plottable CSG type routes through `plottype`, which is what keeps `plot(obj)`
        # working at the doc call sites (D4).
        for obj in (box, box - cone, CSG.surfaces(cone)[3], CSG.Ellipse{Float32}(r = Float32(1)),
                    CartesianPoint{Float32}(1, 2, 3), CylindricalPoint{Float32}(1, Float32(π / 4), 2),
                    [CartesianPoint{Float32}(0, 0, 0)], CartesianVector{Float32}(1, 1, 1))
            @test Makie.plottype(obj) === ext.CSGPlot
        end
    end

    @testset "legend entries" begin
        # Plots labelled the first series of each geometry and let its legend collect it; a Makie
        # recipe produces no entries at all (plan §1.3), so the wrapper builds them explicitly —
        # from `get_label_name`, the same function the Plots recipes now call.
        entries(fap) = MakieSpec._legend_entries(
            only(filter(b -> b isa Makie.Legend, fap.figure.content)))
        @test entries(csgplot(box)) == ["Box"]
        @test entries(csgplot(box - cone)) == ["CSG"]
        @test entries(csgplot(CSG.surfaces(cone)[3])) == ["Cone Mantle"]
        @test entries(csgplot(CSG.Ellipse{Float32}(r = Float32(1)))) == ["Ellipse"]
        # A point's Plots label is the auto-generated `"y1"`, which is not a legend entry on either
        # side, so the figure carries no legend at all.
        @test isempty(filter(b -> b isa Makie.Legend,
                             csgplot(CartesianPoint{Float32}(1, 2, 3)).figure.content))
    end

    @testset "points and vectors" begin
        # The adaptive marker size of the Plots recipe is backend-agnostic and kept, up to the one
        # constant that reconciles the two marker scales.
        single = only_of(csgplot(CartesianPoint{Float32}(1, 2, 3)), :scatter)
        @test single.markersize[][1] ≈ 4 * ext.SSD_MARKERSIZES.point_scale
        many = only_of(csgplot([CartesianPoint{Float32}(cos(φ), sin(φ), 0)
                                for φ in range(0, 2π, length = 17)]), :scatter)
        @test many.markersize[][1] < single.markersize[][1]

        # A cylindrical point is drawn where its Cartesian equivalent is, in display units.
        cyl = only_of(csgplot(CylindricalPoint{Float32}(1, Float32(π / 4), 2)), :scatter)[1][]
        @test only(cyl) ≈ Makie.Point3f(1000cos(π / 4), 1000sin(π / 4), 2000) rtol = 1e-5

        # `plot(::CartesianVector)` has never worked under Plots: its `:vector` series type runs
        # after the positional arguments have been converted to series data, and there is no such
        # conversion for `(CartesianPoint, CartesianVector)`, so the hand-built 13-point arrowhead
        # was unreachable (plan F2). There is therefore no golden spec to gate it against, and this
        # is the assertion instead.
        @test read_spec("csg_vector").note == PlotParityCases.NOTE_VECTOR_BROKEN
        arrow = csgplot(CartesianVector{Float32}(1, 1, 1))
        @test arrow isa Makie.FigureAxisPlot
        @test !isempty(leaves(arrow))
        @test all(k -> k === :meshscatter, keys_of(arrow))   # what `arrows3d!` decomposes into
    end
end

@testset "Detector geometry" begin
    ext = ssd_makie_ext
    det = fixtures().cyl.detector
    leaves(fap) = MakieSpec.leaf_plots(fap isa Makie.FigureAxisPlot ? fap.plot : fap)
    meshes_of(fap) = filter(l -> Makie.plotkey(l) === :mesh, leaves(fap))

    @testset "tier-1 parity gate" begin
        gated = 0
        for c in cases(groups = ["detector"])
            golden = read_spec(c.id)
            isempty(golden.note) && isempty(golden.error) || continue
            args = c.args()
            fig = detectorplot(args...; c.kwargs...)
            spec = makie_case_spec(fig; id = c.id, group = c.group, call = call_string(c, args))
            diffs, _ = case_diffs(golden, spec; ext)
            @test (c.id, diffs) == (c.id, String[])
            gated += 1
        end
        @test gated == 18   # every detector case has a valid Plots reference
    end

    @testset "parts, colours and opacities" begin
        # One list fixes which parts are drawn and in what order, and the legend reads from the same
        # one — so a legend can never name a part the picture does not show.
        @test ext.detector_parts(det) == vcat(collect(det.contacts), collect(det.passives))
        @test ext.detector_parts(det; show_passives = false) == collect(det.contacts)
        @test det.semiconductor in ext.detector_parts(det; show_semiconductor = true)
        # Contacts come first, then the semiconductor, then the passives.
        @test ext.detector_parts(det; show_semiconductor = true)[length(det.contacts) + 1] ===
              det.semiconductor

        # A contact takes its palette colour, the semiconductor grey, a passive silver.
        @test ext.part_color(det.contacts[1]) == ext.contact_color(1)
        @test ext.part_color(det.contacts[2]) == ext.contact_color(2)
        @test ext.part_color(det.semiconductor) == ext.SSD_COLORS.semiconductor
        @test ext.part_color(det.passives[1]) == ext.SSD_COLORS.passive

        # The semiconductor is drawn at 0.5 where contacts and passives are at 0.2: the Plots
        # `Semiconductor` recipe is the one that sets no `fillalpha`, so it keeps the `Mesh` recipe's
        # own default. The goldens record the difference, and it is why a semiconductor reads as the
        # body with the contacts as faint shells on it.
        @test ext.part_alpha(det.semiconductor) == 0.5
        @test ext.part_alpha(det.contacts[1]) == ext.part_alpha(det.passives[1]) == 0.2
        @test only(meshes_of(detectorplot(det.contacts[1]))).alpha[] ≈ 0.2
        @test first(meshes_of(detectorplot(det.semiconductor))).alpha[] ≈ 0.5
        @test only(meshes_of(detectorplot(det.contacts[2]))).color[] == ext.contact_color(2)

        # A detector part draws its wireframe in its own colour, where a bare CSG geometry draws it
        # black — that is what tells two contacts apart in a `:wireframe` plot.
        wires = filter(l -> Makie.plotkey(l) === :linesegments,
                       leaves(detectorplot(det.contacts[2]; seriestype = :wireframe)))
        @test all(w -> w.color[] == ext.contact_color(2), wires)
        @test only(filter(l -> Makie.plotkey(l) === :linesegments,
                          leaves(csgplot(det.contacts[2].geometry;
                                         seriestype = :wireframe)))).color[] ==
              Makie.to_color(:black)
    end

    @testset "one sampling scale for the whole detector" begin
        # The marker representations divide `CSG_scale` by `n_samples`, so sampling each part
        # against its own extent would draw the small contact ten times finer than the crystal and
        # it would come out a solid blob. Both backends take the maximum over all parts, which means
        # a part drawn inside the detector gets *fewer* points than the same part drawn alone.
        n_points(fap) = sum(length(l[1][]) for l in leaves(fap) if Makie.plotkey(l) === :scatter)
        alone = n_points(detectorplot(det.contacts[1]; seriestype = :samplesurface))
        scales = [SolidStateDetectors.ConstructiveSolidGeometry.get_scale(p.geometry)
                  for p in ext.detector_parts(det)]
        @test ext.part_color(det.contacts[1]) !== nothing        # contact 1 is not the largest part
        @test SolidStateDetectors.ConstructiveSolidGeometry.get_scale(det.contacts[1].geometry) <
              maximum(scales)

        whole = detectorplot(det; seriestype = :samplesurface)
        contact1_in_detector = length(first(filter(l -> Makie.plotkey(l) === :scatter,
                                                   leaves(whole)))[1][])
        @test contact1_in_detector < alone
    end

    @testset "entry points" begin
        fap = Makie.plot(det)
        @test fap isa Makie.FigureAxisPlot
        @test fap.plot isa ext.DetectorPlot
        @test fap.axis isa Makie.Axis3
        @test (fap.axis.xlabel[], fap.axis.zlabel[]) == ("x / mm", "z / mm")

        fig = Makie.Figure()
        ax = Makie.Axis3(fig[1, 1])
        @test detectorplot!(ax, det) isa ext.DetectorPlot
        @test ax.title[] == ""
        @test Makie.plot!(ax, det) isa ext.DetectorPlot

        for obj in (det, det.contacts[1], det.semiconductor, det.passives[1])
            @test Makie.plottype(obj) === ext.DetectorPlot
        end
        @test Makie.preferred_axis_type(Makie.plot(det, seriestype = :slice, φ = 0).plot) ===
              Makie.Axis
        @test Makie.preferred_axis_type(Makie.plot(det).plot) === Makie.Axis3
    end

    @testset "legend" begin
        entries(fap) = MakieSpec._legend_entries(
            only(filter(b -> b isa Makie.Legend, fap.figure.content)))
        @test entries(detectorplot(det)) == read_spec("detector_cyl_csg").axes[1].legend_entries
        @test entries(detectorplot(det; show_semiconductor = true)) ==
              read_spec("detector_cyl_csg_semiconductor").axes[1].legend_entries
        @test entries(detectorplot(det; show_passives = false)) ==
              read_spec("detector_cyl_csg_nopassives").axes[1].legend_entries
        @test entries(detectorplot(det.contacts[1])) ==
              read_spec("detector_cyl_contact1").axes[1].legend_entries

        # The swatch follows the representation: a patch for a mesh, a line for a wireframe, a
        # marker for the two marker-based ones.
        _, labels = ext.legend_entries(det)
        @test labels == read_spec("detector_cyl_csg").axes[1].legend_entries
        @test all(e -> e isa Makie.PolyElement, ext.legend_entries(det)[1])
        @test all(e -> e isa Makie.LineElement,
                  ext.legend_entries(det; seriestype = :wireframe)[1])
        @test all(e -> e isa Makie.MarkerElement,
                  ext.legend_entries(det; seriestype = :samplesurface)[1])
    end
end

@testset "Drift paths and charge clouds" begin
    ext = ssd_makie_ext
    CCM = SolidStateDetectors
    leaves(fap) = MakieSpec.leaf_plots(fap isa Makie.FigureAxisPlot ? fap.plot : fap)
    keys_of(fap) = [Makie.plotkey(l) for l in leaves(fap)]
    of_key(fap, k) = filter(l -> Makie.plotkey(l) === k, leaves(fap))
    legend_of(fap) = filter(b -> b isa Makie.Legend, fap.figure.content)
    entries(fap) = MakieSpec._legend_entries(only(legend_of(fap)))

    @testset "tier-1 parity gate" begin
        gated = 0
        for c in cases(groups = ["driftpaths", "chargeclouds"])
            golden = read_spec(c.id)
            isempty(golden.note) && isempty(golden.error) || continue
            args = c.args()
            entry = c.group == "driftpaths" ? driftpathplot : chargecloudplot
            fig = entry(args...; c.kwargs...)
            spec = makie_case_spec(fig; id = c.id, group = c.group, call = call_string(c, args))
            diffs, _ = case_diffs(golden, spec; ext)
            @test (c.id, diffs) == (c.id, String[])
            gated += 1
        end
        @test gated == 13   # 3 drift-path cases and 10 charge-cloud ones
    end

    @testset "drift paths" begin
        evt = fixtures(with_event = true).evt
        paths = evt.drift_paths

        # One `lines!` and one `scatter!` per species, whatever the number of paths: the Plots recipe
        # emits four series *per path*, which for a thousand carriers is four thousand series.
        @test keys_of(driftpathplot(paths)) == [:lines, :scatter, :lines, :scatter]
        @test keys_of(driftpathplot(paths[1])) == keys_of(driftpathplot(paths))

        # The paths are joined by a non-finite vertex so Makie lifts the pen between carriers rather
        # than drawing a line from where one was collected to where the next one started.
        electrons = first(of_key(driftpathplot(paths), :lines))[1][]
        n_breaks = count(p -> !all(isfinite, p), electrons)
        @test n_breaks == length(paths) - 1
        # …and one end marker per carrier, at the last point of its trajectory.
        @test length(first(of_key(driftpathplot(paths), :scatter))[1][]) == length(paths)

        e_line, h_line = of_key(driftpathplot(paths), :lines)
        @test e_line.color[] == ext.SSD_COLORS.electron
        @test h_line.color[] == ext.SSD_COLORS.hole

        # `showlabel` now only decides whether the legend is drawn: a Makie legend is built from
        # explicit elements, so there is one entry per species by construction.
        @test entries(driftpathplot(paths)) == ["e-", "h+"]
        @test entries(driftpathplot(paths)) == read_spec("driftpaths_single").axes[1].legend_entries
        @test isempty(legend_of(driftpathplot(paths; showlabel = false)))

        fap = Makie.plot(paths)
        @test fap.plot isa ext.DriftPathPlot
        @test fap.axis isa Makie.Axis3
        @test (fap.axis.xlabel[], fap.axis.zlabel[]) == ("x / mm", "z / mm")
        fig = Makie.Figure(); ax = Makie.Axis3(fig[1, 1])
        @test driftpathplot!(ax, paths) isa ext.DriftPathPlot
        @test ax.title[] == ""
    end

    @testset "charge clouds" begin
        center = CartesianPoint{Float32}(0, 0, 0)
        ico = NBodyChargeCloud(center, Float32(1) * u"eV";
                               number_of_shells = 2, shell_structure = CCM.Icosahedron)
        shell = CCM.Icosahedron(center, Float32(1))

        # A bare shell is one group; a cloud is its center plus one group per shell.
        @test length(ext._cloud_groups(shell, 12)) == 1
        groups = ext._cloud_groups(ico, 20)
        @test length(groups) == 3                       # center + 2 shells
        @test length(groups[1].points) == 1
        @test all(g -> length(g.points) == 12, groups[2:3])   # icosahedron: 12 vertices
        # Shells shrink by a factor e from the center outwards, so the outer ones read as a halo.
        @test groups[2].markersize ≈ 20 / exp(1)
        @test groups[3].markersize ≈ 20 / exp(2)

        # The whole cloud is two plot objects: one scatter with a size per point, one line with the
        # shell edges. The Plots recipe emits 2 + 2 x 10 series for this cloud.
        @test keys_of(chargecloudplot(ico)) == [:scatter, :lines]
        sizes = only(of_key(chargecloudplot(ico), :scatter)).markersize[]
        @test length(unique(sizes)) == 3                 # one size per group
        @test keys_of(chargecloudplot(ico; connect = false)) == [:scatter]

        # Both backends read the same connection table, so they cannot wire a solid up differently.
        @test SolidStateDetectors.charge_cloud_connections(shell) ==
              SolidStateDetectors.PLATONIC_SOLID_CONNECTIONS[12]
        @test isempty(SolidStateDetectors.charge_cloud_connections(CCM.PointCharge(center, Float32(1))))
        @test keys_of(chargecloudplot(CCM.PointCharge(center, Float32(1)))) == [:scatter]

        # The `N`-charges-per-shell constructor distributes points over a sphere with no structure,
        # so there is nothing to connect.
        plain = NBodyChargeCloud(center, Float32(1) * u"eV", 40; number_of_shells = 2)
        @test keys_of(chargecloudplot(plain)) == [:scatter]

        @test entries(chargecloudplot(ico)) ==
              read_spec("chargecloud_nbody_shells").axes[1].legend_entries
        @test entries(chargecloudplot(shell)) == ["Icosahedron"]
    end

    @testset "particle sources" begin
        pos = CartesianPoint(0.04, 0.0, 0.05)
        # An isotropic source has no direction to draw, so it is its position alone.
        iso = SolidStateDetectors.MonoenergeticSource("gamma", 1u"MeV", pos)
        @test keys_of(chargecloudplot(iso)) == [:scatter]
        @test entries(chargecloudplot(iso)) ==
              read_spec("chargecloud_monoenergeticsource").axes[1].legend_entries

        # A perfectly directed source gets a ray, a partly directed one an emission cone.
        ray = SolidStateDetectors.MonoenergeticSource("gamma", 1u"MeV", pos,
                                                      CartesianVector(0, 0, 1), 0u"°")
        @test :lines in keys_of(chargecloudplot(ray))
        cone = SolidStateDetectors.MonoenergeticSource("gamma", 1u"MeV", pos,
                                                       CartesianVector(0, 0, 1), 30u"°")
        @test :mesh in keys_of(chargecloudplot(cone))
        # …and one wider than a hemisphere is drawn as its position alone: a cone that big would
        # hide the geometry around it rather than describe the source.
        wide = SolidStateDetectors.MonoenergeticSource("gamma", 1u"MeV", pos,
                                                       CartesianVector(0, 0, 1), 150u"°")
        @test keys_of(chargecloudplot(wide)) == [:scatter]

        @test Makie.plot(iso).plot isa ext.ChargeCloudPlot
        @test Makie.plottype(iso) === ext.ChargeCloudPlot
    end
end


@testset "Waveforms and grids" begin
    ext = ssd_makie_ext
    leaves(fap) = MakieSpec.leaf_plots(fap isa Makie.FigureAxisPlot ? fap.plot : fap)
    keys_of(fap) = [Makie.plotkey(l) for l in leaves(fap)]
    legend_of(fig) = filter(b -> b isa Makie.Legend, fig.content)
    entries(fap) = MakieSpec._legend_entries(only(legend_of(fap.figure)))
    axes_of(fig) = [b for b in fig.content if b isa Makie.Axis]

    @testset "tier-1 parity gate" begin
        gated = 0
        for c in cases(groups = ["waveforms", "grids"])
            golden = read_spec(c.id)
            isempty(golden.note) && isempty(golden.error) || continue
            args = c.args()
            # `plot(evt, sim, id)` is not an entry point: the electron/hole contribution is a named
            # function on both sides (D9). Everything else goes through the generic `Makie.plot`.
            fig = c.kind === :ehcontribution ? ehcontributionplot(args...; c.kwargs...) :
                                               Makie.plot(args...; c.kwargs...)
            spec = makie_case_spec(fig; id = c.id, group = c.group, call = call_string(c, args))
            diffs, _ = case_diffs(golden, spec; ext)
            @test (c.id, diffs) == (c.id, String[])
            gated += 1
        end
        @test gated == 7   # 5 waveform cases and 2 grid ones
    end

    @testset "waveforms" begin
        evt = fixtures(with_event = true).evt

        # One `lines!` per waveform, and the labels RadiationDetectorSignals' Plots recipe sets.
        @test keys_of(waveformplot(evt.waveforms[1])) == [:lines]
        @test keys_of(waveformplot(evt.waveforms)) == [:lines, :lines]
        fap = waveformplot(evt.waveforms[1])
        @test (fap.axis.xlabel[], fap.axis.ylabel[]) == ("t / ns", "signal / e")
        @test isempty(legend_of(fap.figure))

        # Colour by position, not by order of appearance: a contact whose waveform is `missing` must
        # not recolour the ones after it, which is what the Plots recipe's `skipmissing` does.
        with_gap = Union{Missing, typeof(evt.waveforms[1])}[missing, evt.waveforms[2]]
        drawn = only(leaves(waveformplot(with_gap)))
        @test drawn.color[] == ext.contact_color(2)

        @test Makie.plot(evt.waveforms[1]).plot isa ext.WaveformPlot
        @test Makie.plottype(evt.waveforms) === ext.WaveformPlot
        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1])
        @test waveformplot!(ax, evt.waveforms) isa ext.WaveformPlot
        @test ax.ylabel[] == ""
    end

    @testset "electron and hole contribution" begin
        evt, sim = fixtures(with_event = true).evt, fixtures().cyl

        fap = ehcontributionplot(evt, sim, 1)
        @test keys_of(fap) == [:lines, :lines, :lines]
        @test entries(fap) == read_spec("waveform_ehcontribution").axes[1].legend_entries
        e_line, h_line, sum_line = leaves(fap)
        @test e_line.color[] == ext.SSD_COLORS.electron
        @test h_line.color[] == ext.SSD_COLORS.hole
        @test sum_line.color[] == ext.contact_color(1)
        @test sum_line.alpha[] == 0.5

        # F4: the Plots recipe writes the *boolean* `false` into its y guide, because `-->` binds
        # tighter than `?:`. Naming the quantity is the accepted deviation (plan §2.7).
        @test read_spec("waveform_ehcontribution").axes[1].yguide == "false"
        @test fap.axis.ylabel[] == "Charge / e"

        # All three curves are extended to one length, so they can be read against each other.
        n = length(evt.waveforms[1].signal)
        @test all(l -> length(l[1][]) == n, leaves(fap))
        @test all(l -> length(l[1][]) == 300, leaves(ehcontributionplot(evt, sim, 1; n_samples = 300)))

        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1])
        @test ehcontributionplot!(ax, evt, sim, 1) isa ext.EHContributionPlot
        @test isempty(legend_of(fig))
    end

    @testset "grids" begin
        cyl = fixtures().cyl.electric_potential.grid
        car = fixtures().car.electric_potential.grid

        # One panel per axis with more than one tick: the cylindrical fixture is 2-D in φ.
        @test ext.grid_plot_dims(cyl) == [1, 3]
        @test ext.grid_plot_dims(car) == [1, 2, 3]
        fig = gridplot(cyl)
        @test fig isa Makie.Figure          # a GridLayout has no single axis to hand back
        panels = axes_of(fig)
        @test length(panels) == 2
        @test [a.title[] for a in panels] ==
              [SolidStateDetectors.grid_panel_title(d) for d in (1, 3)]
        @test [a.xlabel[] for a in panels] == ["r / mm", "z / mm"]
        @test all(a -> a.ylabel[] == ext.GRID_COUNT_NAME, panels)
        @test length(axes_of(gridplot(car))) == 3

        # The histogram is taken on the ticks as the grid holds them and only the edges are
        # rescaled; binning the rescaled ticks moves one across a bin edge and loses a count.
        edges, counts = ext.grid_point_histogram(cyl, 1)
        @test length(edges) == length(counts) + 1
        @test sum(counts) == length(cyl[1].ticks)
        @test maximum(counts) == read_spec("grid_cyl").series[1].y.extrema[2]
        @test extrema(edges) == Tuple(read_spec("grid_cyl").series[1].x.extrema .* 1000)

        # A count axis gets whole-number ticks: Makie's automatic ones put 2.5 on an axis whose
        # values are grid points.
        @test ext.count_ticks(7) == [0, 2, 4, 6]
        @test ext.count_ticks(4) == [0, 1, 2, 3, 4]
        @test all(t -> t == round(t), panels[1].yticks[])

        # The outline closes at zero on both ends, which is what makes it read as a histogram.
        pts = ext.step_outline(edges, counts)
        @test (pts[1][2], pts[end][2]) == (0, 0)
        @test length(pts) == 2 * length(counts) + 2

        @test Makie.plottype(cyl) === ext.GridPlot
        fig2 = Makie.Figure(); ax = Makie.Axis(fig2[1, 1])
        @test gridplot!(ax, cyl; dim = 3) isa ext.GridPlot
        @test ax.title[] == ""
    end
end


@testset "Field lines" begin
    ext = ssd_makie_ext
    leaves(fap) = MakieSpec.leaf_plots(fap isa Makie.FigureAxisPlot ? fap.plot : fap)
    keys_of(fap) = [Makie.plotkey(l) for l in leaves(fap)]

    @testset "tier-1 parity gate" begin
        gated = 0
        for c in cases(groups = ["fieldlines"])
            golden = read_spec(c.id)
            isempty(golden.note) && isempty(golden.error) || continue
            fig = fieldlinesplot(c.args()...; c.kwargs...)
            spec = makie_case_spec(fig; id = c.id, group = c.group, call = call_string(c))
            diffs, _ = case_diffs(golden, spec; ext)
            @test (c.id, diffs) == (c.id, String[])
            gated += 1
        end
        @test gated == 11   # 12 cases minus the radial one Plots' `r` alias makes unusable
    end

    @testset "lines" begin
        sim = fixtures().cyl

        # Every line is one `lines!` with a non-finite vertex between them: the Plots recipe emits
        # one series per line, which for this detector is 96.
        fap = fieldlinesplot(sim; φ = 0, sampling = 8u"mm")
        @test keys_of(fap) == [:lines]
        pts = only(leaves(fap))[1][]
        n_breaks = count(p -> !all(isfinite, p), pts)
        @test n_breaks == 95                      # 96 lines, 95 separators
        @test length(read_spec("fieldlines_cyl_phi").series) == 96

        # `spacing` thins the lines by index, so it must halve them rather than change their shape.
        thinned = only(leaves(fieldlinesplot(sim; φ = 0, sampling = 8u"mm", spacing = 2)))[1][]
        @test count(p -> !all(isfinite, p), thinned) == 47

        # White when drawn into an existing axis, which is how they go over a field-magnitude plot.
        # The form that builds its own figure has nothing underneath and draws them dark: white
        # lines on a white axis are an empty picture, which is what the Plots recipe produced.
        over = Makie.Figure()
        over_ax = Makie.Axis(over[1, 1])
        @test only(leaves(fieldlinesplot!(over_ax, sim; φ = 0, sampling = 8u"mm"))).color[] ==
              ext.SSD_COLORS.fieldline
        @test only(leaves(fap)).color[] == ext.SSD_COLORS.fieldline_alone
        @test only(leaves(fieldlinesplot(sim; φ = 0, sampling = 8u"mm", color = :red))).color[] ==
              Makie.to_color(:red)
    end

    @testset "cross sections" begin
        sim = fixtures().cyl
        # A cylindrical z cut is polar, and a PolarAxis takes its angle in radians — the one place
        # the display-unit convention does not apply (D7).
        polar = fieldlinesplot(sim; z = 0.04, sampling = 8u"mm")
        @test polar.axis isa Makie.PolarAxis
        θ = [p[1] for p in only(leaves(polar))[1][] if all(isfinite, p)]
        @test maximum(abs, θ) <= 2π + 1e-6

        # `full_det` annotates the title only: the sample lines already run from -r to +r, so a
        # field line reaches both halves either way.
        plain = fieldlinesplot(sim; φ = 0, sampling = 8u"mm")
        full = fieldlinesplot(sim; φ = 0, sampling = 8u"mm", full_det = true)
        @test occursin("on left side", full.axis.title[])
        @test !occursin("on left side", plain.axis.title[])
        @test length(only(leaves(full))[1][]) == length(only(leaves(plain))[1][])
    end

    @testset "a cross section no contact meets draws nothing" begin
        # The CGD contacts are slabs perpendicular to x, so no contact surface intersects an x cross
        # section and no charge is ever spawned. Drawing nothing without erroring is the parity
        # property here, and the harness must report no series rather than mistaking the recipe
        # itself for a primitive.
        fap = fieldlinesplot(fixtures().car; x = 0, sampling = 1u"mm")
        @test isempty(leaves(fap))
        spec = makie_case_spec(fap; id = "fieldlines_car_x", group = "fieldlines", call = "")
        @test isempty(spec.series)
        @test isempty(read_spec("fieldlines_car_x").series)
    end

    @testset "F5: a cross section cannot have an axis in the dimension it cut" begin
        # `plot_electric_fieldlines` sets `yguide --> "z"` whatever the cut, so a Cartesian z cross
        # section labels its y data "z". The Makie labels follow the data.
        golden = read_spec("fieldlines_car_z")
        @test (golden.axes[1].xguide, golden.axes[1].yguide) == ("x", "z")
        fap = fieldlinesplot(fixtures().car; z = 0, sampling = 1u"mm")
        @test (fap.axis.xlabel[], fap.axis.ylabel[]) == ("x / mm", "y / mm")
        @test _cut_dimension(golden.axes[1].title) == "z"
    end

    @testset "the shared computation is one algorithm" begin
        sim = fixtures().cyl
        # Both backends read the spawn positions, the drift and the projection from `plot_utils.jl`,
        # so they cannot trace different lines. The ordering of the spawn points is what `spacing`
        # thins, so it is part of the contract.
        pos = SolidStateDetectors.fieldline_spawn_positions(sim, :φ, 0.0f0; sampling = 8u"mm")
        @test length(pos) == 48
        # `((spacing - 1) + i) % spacing == 0` keeps the *odd* indices at `spacing = 2`.
        every_second = SolidStateDetectors.fieldline_spawn_positions(sim, :φ, 0.0f0;
                                                                    sampling = 8u"mm", spacing = 2)
        @test every_second == pos[1:2:end]
    end
end
