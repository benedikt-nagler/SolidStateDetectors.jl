# Makie recipes and figure-level wrappers for waveforms: one `RDWaveform`, the vector of them an
# `Event` carries, and the electron/hole decomposition of one contact's signal.
# 5 of the ~142 parity cases.
#
# The single-waveform picture is not this package's under Plots at all: `plot(::RDWaveform)` is a
# recipe of RadiationDetectorSignals.jl, which has no Makie half. So the guides it sets ("t",
# "signal") are reproduced here rather than shared through `plot_utils.jl` the way the labels this
# package owns are. The `waveform_single` golden spec is what pins them.
#
# One deliberate departure from the Plots recipes: the vector form colours each waveform by its
# *position*, i.e. by contact id. The Plots `Vector{Union{Missing, RDWaveform}}` recipe re-emits
# `collect(skipmissing(wvs))`, so a detector whose first contact has no waveform silently recolours
# every other contact. Indexing by position keeps a contact's curve the colour that contact is drawn
# in everywhere else.

# The two guides RadiationDetectorSignals' Plots recipe sets.
const WAVEFORM_TIME_NAME = "t"
const WAVEFORM_SIGNAL_NAME = "signal"

"""
    waveformplot(waveforms; kwargs...)
    waveformplot!(ax, waveforms; kwargs...)

Plot one `RDWaveform`, or every waveform of an `Event`, against time.

The non-mutating form returns a `FigureAxisPlot` carrying the axis labels. The mutating form draws
into an existing axis and sets none. `plot(waveforms)` routes here.

A `missing` entry is a contact whose signal was not calculated. It is skipped, and the contacts that
do have a waveform keep their own colour rather than closing the gap.

## Keywords
* `color`: colour of the curves. `automatic` gives waveform `i` the colour of contact `i`; a vector
  is cycled over the waveforms.
* `linewidth`, `linestyle`, `alpha`: the curves.
"""
@recipe WaveformPlot (waveforms,) begin
    "Colour of the curves; `automatic` gives waveform `i` the colour of contact `i`."
    color = Makie.automatic
    "Width of the curves."
    linewidth = Makie.automatic
    "Style of the curves."
    linestyle = :solid
    "Opacity of the curves."
    alpha = 1.0
end

const WaveformLike = Union{RDWaveform, AbstractVector{<:Union{Missing, RDWaveform}}}

for WT in (:RDWaveform, :(AbstractVector{<:Union{Missing, RDWaveform}}))
    @eval Makie.plottype(::$WT) = WaveformPlot
    @eval Makie.convert_arguments(::Type{<:WaveformPlot}, w::$WT) = (w,)
end


## ----------------------------------------------------------------------------------------------
## The recipe
## ----------------------------------------------------------------------------------------------

_waveforms(w::RDWaveform) = [w]
_waveforms(ws::AbstractVector) = ws

# The samples as bare numbers in display units, which is how every recipe here passes data (D7).
# `display_unit` leaves the two units a waveform actually carries alone: nanoseconds are already the
# display unit for a time, and a charge or energy calibration is at a sensible scale by construction.
# So this is a no-op on the usual inputs and a rescale on a waveform timed in seconds.
waveform_time(wf) = display_values(Unitful.ustrip.(wf.time), Unitful.unit(eltype(wf.time)))
waveform_signal(wf) = display_values(Unitful.ustrip.(wf.signal), Unitful.unit(eltype(wf.signal)))

_waveform_color(::Makie.Automatic, i::Integer) = contact_color(i)
_waveform_color(c::AbstractVector, i::Integer) = c[mod1(i, length(c))]
_waveform_color(c, ::Integer) = c

function Makie.plot!(p::WaveformPlot)
    lw = color_or_default(p.linewidth[], SSD_LINEWIDTHS.waveform)
    for (i, wf) in enumerate(_waveforms(p[1][]))
        ismissing(wf) && continue
        Makie.lines!(p, waveform_time(wf), waveform_signal(wf);
                     color = _waveform_color(p.color[], i), linewidth = lw,
                     linestyle = p.linestyle[], alpha = p.alpha[])
    end
    return p
end


## ----------------------------------------------------------------------------------------------
## Figure-level wrapper
## ----------------------------------------------------------------------------------------------

_first_waveform(w::RDWaveform) = w
function _first_waveform(ws::AbstractVector)
    for w in ws
        ismissing(w) || return w
    end
    nothing
end

# The axis a waveform wants: time against signal, each labelled with the unit its samples carry. A
# vector of waveforms is described by its first present entry, all of them being sampled and
# calibrated the same way.
function _waveform_axis_setup(waveforms, yname::AbstractString = WAVEFORM_SIGNAL_NAME)
    wf = _first_waveform(waveforms)
    (; xname = WAVEFORM_TIME_NAME, yname,
       xunit = wf === nothing ? nothing : Unitful.unit(eltype(wf.time)),
       yunit = wf === nothing ? nothing : Unitful.unit(eltype(wf.signal)))
end

function waveformplot(gp::Makie.GridPosition, waveforms::WaveformLike;
                      axis = (;), legend = (;), labels = nothing, kwargs...)
    setup = _waveform_axis_setup(waveforms)
    if labels === nothing || legend === false
        ax = ssd_axis(gp; setup..., axis...)
        return Makie.AxisPlot(ax, waveformplot!(ax, waveforms; kwargs...))
    end
    # A labelled plot needs room for the legend beside the axis, so it gets a nested layout; an
    # unlabelled one is a bare axis, which is what keeps `waveform_single` at parity (its Plots
    # legend reads "y1" and carries no information).
    gl = panel_layout(gp)
    ax = ssd_axis(gl[1, 1]; setup..., axis...)
    plt = waveformplot!(ax, waveforms; kwargs...)
    add_legend!(gl, waveform_legend_entries(labels;
                                            color = get(kwargs, :color, Makie.automatic),
                                            linewidth = color_or_default(
                                                get(kwargs, :linewidth, Makie.automatic),
                                                SSD_LINEWIDTHS.waveform))...;
                corner = best_legend_corner(plt), legend...)
    Makie.AxisPlot(ax, plt)
end

function waveformplot(waveforms::WaveformLike; figure = (;), kwargs...)
    fig = Makie.Figure(; figure...)
    ap = waveformplot(fig[1, 1], waveforms; kwargs...)
    Makie.FigureAxisPlot(fig, ap.axis, ap.plot)
end

Makie.plot(waveforms::WaveformLike; kwargs...) = waveformplot(waveforms; kwargs...)


## ----------------------------------------------------------------------------------------------
## Electron and hole contribution
## ----------------------------------------------------------------------------------------------

"""
    ehcontributionplot(event, sim, contact_id; kwargs...)
    ehcontributionplot!(ax, event, sim, contact_id; kwargs...)

Plot the waveform of the contact with `contact_id` together with the electron and hole contributions
to it.

The charges of `event` must already have been drifted in `sim` with
[`drift_charges!`](@ref). The non-mutating form returns a `FigureAxisPlot` carrying the axis labels
and the three-entry legend; the mutating form sets neither.

This is the Makie counterpart of [`plot_electron_and_hole_contribution`](@ref), under a different
name so that both backends stay usable in one session.

## Keywords
* `n_samples`: number of samples the curves are drawn with. Defaults to the length of the contact's
  own waveform; shorter truncates, longer extends with the final value.
* `electron_color`, `hole_color`, `sum_color`: colour of each curve.
* `linewidth`: the curves.
* `sum_alpha`: opacity of the summed waveform, drawn over the two contributions.

See also [`get_electron_and_hole_contribution`](@ref).
"""
@recipe EHContributionPlot (event, simulation, contact_id) begin
    "Number of samples the curves are drawn with; `automatic` is the length of the contact's waveform."
    n_samples = Makie.automatic
    "Colour of the electron contribution."
    electron_color = Makie.automatic
    "Colour of the hole contribution."
    hole_color = Makie.automatic
    "Colour of the summed waveform; `automatic` is the contact's own colour."
    sum_color = Makie.automatic
    "Width of the three curves."
    linewidth = Makie.automatic
    "Opacity of the summed waveform, which is drawn over the two contributions."
    sum_alpha = 0.5
end

Makie.convert_arguments(::Type{<:EHContributionPlot}, evt::Event, sim::Simulation, id::Integer) =
    (evt, sim, Int(id))

# The three curves, all extended to the same length so they can be read against each other. The
# electron contribution ends when the last electron is collected, which is well before the hole
# contribution does.
function _eh_waveforms(evt::Event, sim::Simulation, contact_id::Int; n_samples = Makie.automatic)
    wf = evt.waveforms[contact_id]
    signal_unit = Unitful.unit(eltype(wf.signal))
    n = n_samples === Makie.automatic ? length(wf.signal) : Int(n_samples)
    c = get_electron_and_hole_contribution(evt, sim, contact_id; signal_unit)
    (electron = add_baseline_and_extend_tail(c.electron_contribution, 0, n),
     hole     = add_baseline_and_extend_tail(c.hole_contribution, 0, n),
     sum      = add_baseline_and_extend_tail(wf, 0, n))
end

# The colours of the three curves, resolved once so that `ehcontributionplot` builds a legend
# describing exactly what was drawn.
_eh_colors(contact_id::Integer; electron_color = Makie.automatic, hole_color = Makie.automatic,
           sum_color = Makie.automatic, kwargs...) =
    (electron = color_or_default(electron_color, SSD_COLORS.electron),
     hole     = color_or_default(hole_color, SSD_COLORS.hole),
     sum      = color_or_default(sum_color, contact_color(contact_id)))

function Makie.plot!(p::EHContributionPlot)
    contact_id = p[3][]
    wfs = _eh_waveforms(p[1][], p[2][], contact_id; n_samples = p.n_samples[])
    colors = _eh_colors(contact_id; electron_color = p.electron_color[],
                        hole_color = p.hole_color[], sum_color = p.sum_color[])
    lw = color_or_default(p.linewidth[], SSD_LINEWIDTHS.contribution)
    for (part, alpha) in ((:electron, 1.0), (:hole, 1.0), (:sum, p.sum_alpha[]))
        wf = getfield(wfs, part)
        Makie.lines!(p, waveform_time(wf), waveform_signal(wf);
                     color = getfield(colors, part), linewidth = lw, alpha)
    end
    return p
end

function ehcontributionplot(gp::Makie.GridPosition, evt::Event, sim::Simulation,
                            contact_id::Integer; axis = (;), legend = (;), kwargs...)
    wf = evt.waveforms[contact_id]
    # Naming the quantity is a deviation from the golden spec, not an oversight. The Plots recipe
    # writes the *boolean* `false` into its y guide, because `-->` binds tighter than `?:` (plan F4).
    yname = waveform_signal_label(Unitful.unit(eltype(wf.signal)))
    gl = panel_layout(gp)
    ax = ssd_axis(gl[1, 1]; _waveform_axis_setup(wf, yname)..., axis...)
    plt = ehcontributionplot!(ax, evt, sim, contact_id; kwargs...)
    lw = color_or_default(get(kwargs, :linewidth, Makie.automatic), SSD_LINEWIDTHS.contribution)
    legend === false ||
        add_legend!(gl, eh_contribution_legend_entries(_eh_colors(contact_id; kwargs...);
                                                       linewidth = lw)...;
                    corner = best_legend_corner(plt), legend...)
    Makie.AxisPlot(ax, plt)
end

function ehcontributionplot(evt::Event, sim::Simulation, contact_id::Integer;
                            figure = (;), kwargs...)
    fig = Makie.Figure(; figure...)
    ap = ehcontributionplot(fig[1, 1], evt, sim, contact_id; kwargs...)
    Makie.FigureAxisPlot(fig, ap.axis, ap.plot)
end
