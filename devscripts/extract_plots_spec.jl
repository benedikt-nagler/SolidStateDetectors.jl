#!/usr/bin/env julia
#
# Tier 0 of the Plots.jl -> Makie.jl parity harness.
#
# Walks the case matrix (`test/parity/cases.jl`), plots every case with Plots.jl/GR, introspects the
# resulting plot object *without rendering it*, and writes one JSON per case into
# `test/parity/golden/`. Those JSONs are committed: they are the frozen specification of the
# behaviour the Makie port has to reproduce, and they stay valid after the Plots recipes are deleted
# in step 11.
#
# Usage
#
#     julia --project=devscripts devscripts/extract_plots_spec.jl                  # everything
#     julia --project=devscripts devscripts/extract_plots_spec.jl potentials efield # named groups
#     julia --project=devscripts devscripts/extract_plots_spec.jl --fast            # skip slow cases
#     julia --project=devscripts devscripts/extract_plots_spec.jl --list            # list the matrix
#
# Regenerate after any deliberate change to the Plots recipes, and review the JSON diff — an
# unexplained diff there is a behaviour change, which is exactly what the harness exists to surface.

# GR must not try to open a window; the specs never need a rendered image anyway.
ENV["GKSwstype"] = "100"

using Printf
using Plots
using SolidStateDetectors

include(joinpath(@__DIR__, "..", "test", "parity", "cases.jl"))
include(joinpath(@__DIR__, "plots_spec.jl"))

using .PlotParityCases
using .PlotsSpec
using .PlotsSpec.PlotParitySpec

const GOLDEN_DIR = normpath(joinpath(@__DIR__, "..", "test", "parity", "golden"))

"""
    plot_case(c::Case, args)

Call the Plots entry point this case describes and return the `Plot` object.
"""
function plot_case(c::Case, args)
    if c.kind === :plot
        Plots.plot(args...; c.kwargs...)
    elseif c.kind === :fieldlines
        SolidStateDetectors.plot_electric_fieldlines(args...; c.kwargs...)
    elseif c.kind === :ehcontribution
        SolidStateDetectors.plot_electron_and_hole_contribution(args...; c.kwargs...)
    else
        error("unknown case kind $(c.kind)")
    end
end

function main(argv)
    if "--list" in argv
        for c in cases()
            @printf("%-52s %-12s %-7s %s\n", c.id, c.group, c.slow ? "(slow)" : "",
                    isempty(c.note) ? "" : "(no Plots reference)")
        end
        return 0
    end

    include_slow = !("--fast" in argv)
    groups = filter(a -> !startswith(a, "--"), argv)
    isempty(groups) && (groups = nothing)

    cs = cases(; groups, include_slow)
    @info "extracting golden specs" n_cases = length(cs) groups = something(groups, "all") include_slow

    gr()
    # Building the fixtures runs two full simulations plus one drifted event; do it up front so the
    # per-case timings below measure plotting, not physics.
    t_fix = @elapsed fixtures()
    @info "fixtures built" seconds = round(t_fix, digits = 1)

    mkpath(GOLDEN_DIR)
    failed = String[]
    for (i, c) in enumerate(cs)
        t0 = time()
        spec = try
            args = c.args()
            plt = plot_case(c, args)
            plots_case_spec(plt; c.id, c.group, call = call_string(c, args), c.note)
        catch err
            push!(failed, c.id)
            @warn "case failed" id = c.id exception = (err, catch_backtrace())
            PlotsSpec.error_spec(c.id, c.group, call_string(c), err; c.note)
        end
        write_spec(spec; dir = GOLDEN_DIR)
        @printf("[%3d/%3d] %-52s %-6s %5.1fs  (%d series, %d axes)\n", i, length(cs), c.id,
                !isempty(spec.error) ? "FAIL" : isempty(spec.note) ? "ok" : "ok*", time() - t0,
                length(spec.series), length(spec.axes))
    end

    noted = count(c -> !isempty(c.note), cs)
    # A case that carries a note is one Plots is known not to be able to draw (`csg_vector`), so its
    # failure is the expected outcome and must not make the run look broken.
    noted_ids = Set(c.id for c in cs if !isempty(c.note))
    unexpected = filter(!in(noted_ids), failed)
    @info "done" written = length(cs) ok = length(cs) - length(failed) failed = length(failed) noted dir = GOLDEN_DIR
    noted > 0 && @info "`ok*` cases carry a note: their Plots output is not a valid parity reference \
                        (see the spec's `note` field, PlotParityCases.NOTE_R_ALIAS / NOTE_VECTOR_BROKEN)"
    isempty(failed) || @info "failed cases" ids = failed expected = collect(intersect(failed, noted_ids))
    isempty(unexpected) || @warn "unexpected failures" ids = unexpected
    isempty(unexpected) ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
