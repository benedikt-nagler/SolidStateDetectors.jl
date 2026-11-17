#!/usr/bin/env julia
#
# Tier 2 of the Plots.jl -> Makie.jl parity harness.
#
# Renders cases of the matrix in *both* backends at identical pixel size into
# `gallery/{plots,makie}/<case>.png` and writes `gallery/index.html`, a side-by-side page carrying,
# per case, the call, the tier-1 gate result and any note. This is the human sign-off artifact for
# each porting step; it is deliberately not part of CI (see plan §3, "CI runs tiers 0/1 only").
#
# Usage
#
#     julia --project=devscripts devscripts/plot_gallery.jl                   # everything
#     julia --project=devscripts devscripts/plot_gallery.jl potentials efield # named groups
#     julia --project=devscripts devscripts/plot_gallery.jl --fast            # skip slow cases
#     julia --project=devscripts devscripts/plot_gallery.jl --with-plots      # re-render Plots too
#     julia --project=devscripts devscripts/plot_gallery.jl --plots-only      # skip the Makie half
#     julia --project=devscripts devscripts/plot_gallery.jl --no-gate         # skip tier 1
#     julia --project=devscripts devscripts/plot_gallery.jl --page-only       # rebuild the page only
#
# **A filtered run updates only the cases it rendered and keeps the rest.** The rendered rows are
# persisted in `gallery/rows.json` and the page is rebuilt from that manifest merged with the *full*
# case matrix, so `… plot_gallery.jl efield` refreshes the twelve field panels without throwing away
# the hundred-odd panels rendered by an earlier run — and every case that has never been rendered is
# still listed, which is what makes the page a progress view as well as a review artifact.
#
# The Plots half is frozen. The recipes it was rendered from were removed in step 11 of the
# migration, so a run keeps the Plots panels already on disk and renders only Makie. `--with-plots`
# and `--plots-only` are for a checkout from before step 11.
#
# `gallery/` is git-ignored: these PNGs are review material, not repository content (assumption A2).

ENV["GKSwstype"] = "100"

using Printf
using JSON
using Unitful
using Plots
using CairoMakie
using SolidStateDetectors

include(joinpath(@__DIR__, "..", "test", "parity", "cases.jl"))
include(joinpath(@__DIR__, "..", "test", "parity", "makie_spec.jl"))
using .PlotParityCases
using .MakieSpec
using .MakieSpec.PlotParitySpec

# `gate.jl` is plain code included into a scope that already has `PlotParitySpec`, `Unitful` and the
# extension module — the same way `test/test_plotting.jl` includes it. Running the tier-1 gate here
# as well is what lets the page say *why* two panels differ: a pair that looks wrong while the gate
# passes is a rendering problem (the polar sliver of step 3), a pair that looks right while the gate
# fails is a lost title, colorbar or limit.
const ssd_makie_ext = Base.get_extension(SolidStateDetectors, :SolidStateDetectorsMakieExt)
include(joinpath(@__DIR__, "..", "test", "parity", "gate.jl"))

const GALLERY_DIR = normpath(joinpath(@__DIR__, "..", "gallery"))
const MANIFEST = joinpath(GALLERY_DIR, "rows.json")
const SIZE_PX = (600, 450)

## ----------------------------------------------------------------------------------------------
## Rendering
## ----------------------------------------------------------------------------------------------

function render_plots(c::Case, args, path::AbstractString)
    plt = if c.kind === :plot
        Plots.plot(args...; size = SIZE_PX, c.kwargs...)
    elseif c.kind === :fieldlines
        SolidStateDetectors.plot_electric_fieldlines(args...; size = SIZE_PX, c.kwargs...)
    else
        SolidStateDetectors.plot_electron_and_hole_contribution(args...; size = SIZE_PX, c.kwargs...)
    end
    Plots.savefig(plt, path)
    path
end

# The Makie half goes through the same generic entry point the docs use (`Makie.plot`), so it
# exercises the `plottype`/`convert_arguments` wiring rather than a hand-picked function name. The
# figure is returned as well as saved, because the tier-1 gate introspects the same object that was
# rendered — anything else would compare a picture against a different plot.
function render_makie(c::Case, args, path::AbstractString)
    # `plot_electron_and_hole_contribution` is a Plots `@userplot` with no `plot(...)` form, and its
    # Makie counterpart is a named entry point for the same reason (D9), so this one kind cannot go
    # through the generic path.
    fap = if c.kind === :ehcontribution
        SolidStateDetectors.ehcontributionplot(args...; figure = (; size = SIZE_PX), c.kwargs...)
    elseif c.kind === :fieldlines
        SolidStateDetectors.fieldlinesplot(args...; figure = (; size = SIZE_PX), c.kwargs...)
    else
        Makie.plot(args...; figure = (; size = SIZE_PX), c.kwargs...)
    end
    # `px_per_unit = 1` is what makes the two halves the same pixel size: CairoMakie defaults to 2,
    # and a 2x-scaled Makie panel next to a 1x GR panel is not a comparison a reviewer can trust.
    CairoMakie.save(path, fap isa Makie.FigureAxisPlot ? fap.figure : fap; px_per_unit = 1)
    fap
end

# A type with no Makie recipe yet fails in dispatch (`MethodError`), in `plottype` (a plain
# `ErrorException`, "No recipe for plot with args: …") or in argument conversion (an
# `ArgumentError`). All three mean "not ported yet" and are reported as a status rather than as an
# error, so the page's red boxes stay reserved for genuine breakage.
function _is_not_ported(err)
    err isa MethodError && return true
    msg = lowercase(sprint(showerror, err))
    any(occursin(p, msg) for p in ("no recipe for plot", "cannot convert", "conversion failed"))
end

## ----------------------------------------------------------------------------------------------
## Tier 1, per case
## ----------------------------------------------------------------------------------------------

# Gate verdicts, in the order the page ranks them: a failure is what a reviewer must look at first.
const GATE_FAIL, GATE_PASS, GATE_NOREF, GATE_SKIPPED = "fail", "pass", "noref", "skipped"

"""
    gate_case(c, fap) -> (verdict, diffs, deviations)

Run the tier-1 gate on the figure that was just rendered.

`GATE_NOREF` covers the cases whose Plots output is not a valid reference (plan §3.5, F1/F2) and the
ones with no golden spec at all; `GATE_SKIPPED` means there was no Makie figure to compare.
"""
function gate_case(c::Case, fap)
    fap === nothing && return (GATE_SKIPPED, String[], String[])
    PlotParitySpec.has_spec(c.id) || return (GATE_NOREF, ["no golden spec for this case"], String[])
    golden = read_spec(c.id)
    (isempty(golden.note) && isempty(golden.error)) ||
        return (GATE_NOREF, String[], String[])
    spec = makie_case_spec(fap; id = c.id, group = c.group, call = call_string(c))
    diffs, deviations = case_diffs(golden, spec; ext = ssd_makie_ext)
    (isempty(diffs) ? GATE_PASS : GATE_FAIL, diffs, deviations)
end

## ----------------------------------------------------------------------------------------------
## The manifest
##
## One row per case, so that a filtered run can update its own rows and leave the others standing.
## Rows are plain `Dict{String, Any}` because they round trip through JSON.
## ----------------------------------------------------------------------------------------------

read_manifest() =
    isfile(MANIFEST) ? Dict{String, Any}(JSON.parsefile(MANIFEST)) : Dict{String, Any}()

function write_manifest(rows::AbstractDict)
    open(MANIFEST, "w") do io
        JSON.print(io, rows)
    end
    MANIFEST
end

# A row whose PNG has since been deleted must not keep pointing at it.
_existing(rel) = (rel !== nothing && isfile(joinpath(GALLERY_DIR, rel))) ? rel : nothing

# …and a PNG on disk that no row claims — left by a run that predates the manifest — is still worth
# showing, so it is adopted rather than reported as never rendered.
function _panel(recorded, id::AbstractString, dir::AbstractString)
    found = _existing(recorded)
    found === nothing ? _existing("$dir/$id.png") : found
end

function row_of(c::Case; args, plots, makie, makie_status, errors, gate, diffs, deviations)
    Dict{String, Any}(
        "id" => c.id, "group" => c.group, "call" => call_string(c, args), "note" => c.note,
        "error" => join(errors, "\n"), "plots" => plots, "makie" => makie,
        "makie_status" => makie_status, "gate" => gate, "gate_diffs" => diffs,
        "deviations" => deviations, "rendered" => now_string(),
    )
end

# The placeholder for a case the current run did not touch and no earlier run rendered. It doubles
# as the schema: a row read back from an older manifest is merged over this, so a manifest written
# before a field existed still renders.
pending_row(c::Case) = Dict{String, Any}(
    "id" => c.id, "group" => c.group, "call" => call_string(c), "note" => c.note,
    "error" => "", "plots" => nothing, "makie" => nothing,
    "makie_status" => "not rendered yet", "gate" => GATE_SKIPPED, "gate_diffs" => String[],
    "deviations" => String[], "rendered" => "",
)

## ----------------------------------------------------------------------------------------------
## HTML page
## ----------------------------------------------------------------------------------------------

_esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

const HTML_HEAD = """
<!doctype html><html><head><meta charset="utf-8">
<title>SolidStateDetectors.jl — Plots vs Makie parity gallery</title>
<style>
  body { font-family: system-ui, sans-serif; margin: 2rem; background: #fafafa; color: #111; }
  h1 { font-size: 1.4rem; margin-bottom: .25rem; }
  .meta { color: #666; font-size: .85rem; margin-bottom: 1rem; }
  .bar { position: sticky; top: 0; background: #fafafaee; backdrop-filter: blur(4px);
         padding: .6rem 0; border-bottom: 1px solid #ddd; margin-bottom: 1rem; z-index: 5; }
  .bar button { font: inherit; font-size: .8rem; margin: 0 .3rem .3rem 0; padding: .25rem .6rem;
                border: 1px solid #ccc; border-radius: 999px; background: #fff; cursor: pointer; }
  .bar button.on { background: #111; color: #fff; border-color: #111; }
  .bar .grouplabel { font-size: .75rem; color: #888; text-transform: uppercase;
                     letter-spacing: .05em; margin-right: .4rem; }
  .case { background: #fff; border: 1px solid #ddd; border-radius: 6px;
          margin-bottom: 1.5rem; padding: 1rem; }
  .case h2 { font-size: 1rem; margin: 0 0 .25rem 0; font-family: ui-monospace, monospace; }
  .call { font-family: ui-monospace, monospace; font-size: .8rem; color: #555; margin-bottom: .5rem; }
  .chip { display: inline-block; font-size: .7rem; padding: .1rem .5rem; border-radius: 999px;
          margin-left: .4rem; vertical-align: middle; font-family: system-ui, sans-serif; }
  .chip.pass    { background: #e6f4ea; color: #1e7e34; }
  .chip.fail    { background: #fdecea; color: #b02a37; }
  .chip.noref   { background: #fff8e1; color: #8a6d1b; }
  .chip.skipped { background: #eee;    color: #777; }
  .note { background: #fff8e1; border-left: 3px solid #f0ad4e; padding: .4rem .6rem;
          font-size: .8rem; margin-bottom: .5rem; }
  .err  { background: #fdecea; border-left: 3px solid #d9534f; padding: .4rem .6rem;
          font-size: .8rem; font-family: ui-monospace, monospace; white-space: pre-wrap; }
  .diffs { background: #fdecea; border-left: 3px solid #d9534f; padding: .4rem .6rem;
           font-size: .8rem; margin-bottom: .5rem; }
  .diffs li { font-family: ui-monospace, monospace; }
  .dev { background: #eef3fb; border-left: 3px solid #6c8ebf; padding: .4rem .6rem;
         font-size: .8rem; margin-bottom: .5rem; }
  .pair { display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; align-items: start; }
  .pair figure { margin: 0; }
  .pair figcaption { font-size: .75rem; color: #666; text-transform: uppercase;
                     letter-spacing: .05em; margin-bottom: .25rem; }
  .pair img { width: 100%; border: 1px solid #eee; }
  .missing { display: flex; align-items: center; justify-content: center; height: 200px;
             border: 1px dashed #ccc; color: #999; font-size: .85rem; }
  /* A case whose Plots output is not a valid reference still shows the panel — seeing what Plots
     actually did is the point — but it must not be mistaken for the thing Makie has to match. */
  .noref img { border: 2px dashed #f0ad4e; filter: grayscale(.35) opacity(.8); }
  .noref figcaption { color: #b8860b; }
  .caveat { background: #fff3cd; border: 1px solid #f0ad4e; border-radius: 4px; color: #6b4e00;
            padding: .5rem .7rem; font-size: .85rem; margin-bottom: .4rem; }
  .stale { color: #999; font-size: .7rem; float: right; }
</style></head><body>
"""

# Filtering is client side: 142 cases is small enough to ship as one page, and a reviewer working
# through a porting step wants "show me the gate failures" and "show me this group" without a
# regeneration round trip.
const HTML_SCRIPT = """
<script>
(function () {
  const state = {group: "all", status: "all"};
  function apply() {
    for (const el of document.querySelectorAll(".case")) {
      const okGroup = state.group === "all" || el.dataset.group === state.group;
      const okStatus = state.status === "all" || el.dataset.status === state.status;
      el.style.display = (okGroup && okStatus) ? "" : "none";
    }
    document.getElementById("shown").textContent =
      document.querySelectorAll(".case:not([style*='none'])").length;
  }
  for (const b of document.querySelectorAll(".bar button")) {
    b.addEventListener("click", () => {
      const kind = b.dataset.kind;
      state[kind] = b.dataset.value;
      for (const o of document.querySelectorAll(`.bar button[data-kind="\${kind}"]`))
        o.classList.toggle("on", o === b);
      apply();
    });
  }
  apply();
})();
</script>
"""

# What the Plots panel of a case actually shows when it is not the plot the call asked for. Written
# onto the panel itself: a reviewer comparing two pictures does not read the note above them.
function _reference_caveat(r)
    note = r["note"]
    occursin("rotation", note) && return "Plots ignored <code>r = …</code> (it reads <code>r</code> \
        as <code>rotation</code>) and drew the <b>φ = 0 cut</b> instead. This is a different cross \
        section from the one on the right. The Makie panel is the r cut that was asked for."
    occursin("CartesianVector", note) && return "Plots cannot draw a vector. This panel is an \
        unrelated fallback, not the arrow."
    isempty(note) || return "Not the plot the call asked for. See the note above."
    r["group"] == "fieldlines" && return "Plots draws field lines in white, for use over a field \
        plot. On their own they are white on white, so this panel looks empty. The data is the same."
    ""
end

_chip(gate) = gate == GATE_PASS ? "<span class=\"chip pass\">gate ✓</span>" :
              gate == GATE_FAIL ? "<span class=\"chip fail\">gate ✗</span>" :
              gate == GATE_NOREF ? "<span class=\"chip noref\">no reference</span>" :
              "<span class=\"chip skipped\">not gated</span>"

function _filter_bar(io, groups, counts)
    println(io, "<div class=\"bar\">")
    println(io, "<div><span class=\"grouplabel\">group</span>")
    println(io, "<button data-kind=\"group\" data-value=\"all\" class=\"on\">all</button>")
    for g in groups
        @printf(io, "<button data-kind=\"group\" data-value=\"%s\">%s</button>", _esc(g), _esc(g))
    end
    println(io, "</div><div><span class=\"grouplabel\">status</span>")
    println(io, "<button data-kind=\"status\" data-value=\"all\" class=\"on\">all</button>")
    for (value, label) in ((GATE_FAIL, "gate failures"), (GATE_PASS, "gate passes"),
                           (GATE_NOREF, "no reference"), (GATE_SKIPPED, "not gated"))
        @printf(io, "<button data-kind=\"status\" data-value=\"%s\">%s (%d)</button>",
                value, label, get(counts, value, 0))
    end
    println(io, "</div></div>")
end

function write_html(rows, path::AbstractString)
    groups = unique(r["group"] for r in rows)
    counts = Dict(v => count(r -> r["gate"] == v, rows)
                  for v in (GATE_FAIL, GATE_PASS, GATE_NOREF, GATE_SKIPPED))
    n_both = count(r -> r["plots"] !== nothing && r["makie"] !== nothing, rows)
    n_makie = count(r -> r["makie"] !== nothing, rows)
    n_noref = count(r -> !isempty(r["note"]), rows)

    open(path, "w") do io
        print(io, HTML_HEAD)
        println(io, "<h1>Plots vs Makie parity gallery</h1>")
        @printf(io, "<div class=\"meta\"><span id=\"shown\">%d</span> of %d cases shown · \
                     %d rendered in both backends · %d Makie panels · page written %s</div>\n",
                length(rows), length(rows), n_both, n_makie, now_string())
        n_noref == 0 || @printf(io, "<div class=\"note\">%d of these cases have <b>no valid Plots \
            reference</b>: the call does not reach the Plots recipe at all, so its panel shows a \
            different plot. For the <code>r = …</code> cross sections that is the φ = 0 cut, a \
            different slice through the detector. Each such panel says so. Every other pair is \
            checked by the tier-1 gate to show the same cut, title, limits and data.</div>\n", n_noref)
        _filter_bar(io, groups, counts)

        for r in rows
            @printf(io, "<div class=\"case\" data-group=\"%s\" data-status=\"%s\" data-rendered=\"%s\">\n",
                    _esc(r["group"]), _esc(r["gate"]),
                    (r["plots"] !== nothing && r["makie"] !== nothing) ? "both" : "partial")
            println(io, "<h2>", _esc(r["id"]), " <span style=\"color:#888;font-weight:normal\">(",
                    _esc(r["group"]), ")</span>", _chip(r["gate"]),
                    isempty(r["rendered"]) ? "" :
                        "<span class=\"stale\">rendered $(_esc(r["rendered"]))</span>", "</h2>")
            println(io, "<div class=\"call\">", _esc(r["call"]), "</div>")
            isempty(r["note"]) || println(io, "<div class=\"note\">", _esc(r["note"]), "</div>")
            isempty(r["error"]) || println(io, "<div class=\"err\">", _esc(r["error"]), "</div>")
            if !isempty(r["gate_diffs"])
                println(io, "<div class=\"diffs\"><b>tier-1 gate</b><ul>")
                for d in r["gate_diffs"]
                    println(io, "<li>", _esc(d), "</li>")
                end
                println(io, "</ul></div>")
            end
            if !isempty(r["deviations"])
                println(io, "<div class=\"dev\"><b>accepted deviations</b><ul>")
                for d in r["deviations"]
                    println(io, "<li>", _esc(d), "</li>")
                end
                println(io, "</ul></div>")
            end
            println(io, "<div class=\"pair\">")
            # The Plots panel of a noted case is *not* the reference — say so on the panel itself,
            # not only in the note above it, or a reviewer reads it as one.
            noref = !isempty(r["note"])
            caveat = _reference_caveat(r)
            for (backend, img, status, mark) in
                    (("Plots / GR", r["plots"], "not rendered", noref),
                     ("Makie / CairoMakie", r["makie"], r["makie_status"], false))
                println(io, "<figure", mark ? " class=\"noref\"" : "", ">")
                println(io, "<figcaption>", backend,
                        mark ? " — not a parity reference" : "", "</figcaption>")
                backend == "Plots / GR" && !isempty(caveat) &&
                    println(io, "<div class=\"caveat\">", caveat, "</div>")
                println(io, img === nothing ?
                    "<div class=\"missing\">$(_esc(status))</div>" :
                    "<img src=\"$(_esc(img))\" loading=\"lazy\">")
                println(io, "</figure>")
            end
            println(io, "</div></div>")
        end
        print(io, HTML_SCRIPT)
        println(io, "</body></html>")
    end
    path
end

# `Libc.strftime` rather than `Dates`, to keep this script's environment to the two backends.
now_string() = Base.Libc.strftime("%Y-%m-%d %H:%M", time())

## ----------------------------------------------------------------------------------------------

function render_case(c::Case; plots_only::Bool, gate::Bool, with_plots::Bool = false)
    errs = String[]
    args = nothing
    try
        args = c.args()
    catch err
        push!(errs, "building arguments: " * sprint(showerror, err))
    end

    p_rel = m_rel = nothing
    fap = nothing
    m_status = plots_only ? "skipped (--plots-only)" : "not rendered"
    if args !== nothing
        if with_plots || plots_only
            try
                render_plots(c, args, joinpath(GALLERY_DIR, "plots", "$(c.id).png"))
                p_rel = "plots/$(c.id).png"
            catch err
                push!(errs, "Plots: " * sprint(showerror, err))
            end
        end
        if !plots_only
            try
                fap = render_makie(c, args, joinpath(GALLERY_DIR, "makie", "$(c.id).png"))
                m_rel = "makie/$(c.id).png"
            catch err
                fap = nothing
                if _is_not_ported(err)
                    m_status = "not ported yet"
                else
                    m_status = "render failed"
                    push!(errs, "Makie: " * sprint(showerror, err))
                end
            end
        end
    end

    verdict, diffs, deviations = GATE_SKIPPED, String[], String[]
    if gate && fap !== nothing
        try
            verdict, diffs, deviations = gate_case(c, fap)
        catch err
            verdict = GATE_FAIL
            diffs = ["gate errored: " * sprint(showerror, err)]
        end
    end

    row_of(c; args, plots = p_rel, makie = m_rel, makie_status = m_status, errors = errs,
           gate = verdict, diffs, deviations)
end

function main(argv)
    include_slow = !("--fast" in argv)
    plots_only = "--plots-only" in argv
    with_plots = "--with-plots" in argv
    page_only = "--page-only" in argv
    gate = !("--no-gate" in argv) && !plots_only
    groups = filter(a -> !startswith(a, "--"), argv)
    isempty(groups) && (groups = nothing)

    selected = page_only ? Case[] : cases(; groups, include_slow)
    @info "rendering gallery" n_cases = length(selected) groups = something(groups, "all") plots_only gate

    Plots.gr()
    mkpath(joinpath(GALLERY_DIR, "plots"))
    mkpath(joinpath(GALLERY_DIR, "makie"))

    manifest = read_manifest()
    if !isempty(selected)
        t_fix = @elapsed fixtures()
        @info "fixtures built" seconds = round(t_fix, digits = 1)
    end
    for (i, c) in enumerate(selected)
        row = render_case(c; plots_only, gate, with_plots)
        manifest[c.id] = row
        @printf("[%3d/%3d] %-52s plots:%-3s makie:%-3s gate:%s\n", i, length(selected), c.id,
                row["plots"] === nothing ? "—" : "ok", row["makie"] === nothing ? "—" : "ok",
                row["gate"])
    end
    write_manifest(manifest)

    # The page always lists the *whole* matrix, in matrix order: the rows this run produced, the
    # rows an earlier run left in the manifest, and a placeholder for everything still untouched.
    rows = map(cases()) do c
        row = merge(pending_row(c), Dict{String, Any}(get(manifest, c.id, Dict{String, Any}())))
        row["plots"] = _panel(row["plots"], c.id, "plots")
        row["makie"] = _panel(row["makie"], c.id, "makie")
        row
    end
    path = write_html(rows, joinpath(GALLERY_DIR, "index.html"))
    @info "gallery written" path n_rows = length(rows) n_makie =
        count(r -> r["makie"] !== nothing, rows) n_gate_fail = count(r -> r["gate"] == GATE_FAIL, rows)
    0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
