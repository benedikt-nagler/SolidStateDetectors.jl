# # Inactive layer simulation: True coaxial detector

# ## Set the physical models for the inactive layer (dead layer) simulation with YAML
# - using drift model for the inactive layer
# - using constant lifetime trapping model
# - using mobility-tied diffusion model (by calling the `calculate_mobility` function with the drift model)

using CairoMakie
using Unitful
using SolidStateDetectors
T = Float64;

# the geometry parameters of the model for following display
det_rin = 1u"mm"
det_z = det_r = 10u"mm"
z_draw = det_z/2
pn_r = 8.957282u"mm" # this one was calculated by searching the zero impurity point (displayed in the following section)

sim = Simulation{T}(SSD_examples[:TrueCoaxial])
cfn = SSD_examples[:TrueCoaxial]
print(open(f -> read(f, String), cfn))
# `:wireframe` rather than the default `:csg`: CairoMakie draws a translucent mesh triangle by
# triangle, so a 3-D detector reads as a spiky web (see the warning in the Plotting manual).
detectorplot(sim.detector; seriestype = :wireframe)
#jl save("tutorial_det_dl.pdf", current_figure()) # hide
#md save("tutorial_det_dl.pdf", current_figure()) # hide
#md save("tutorial_det_dl.svg", current_figure()); nothing # hide
#md # [![tutorial_det_dl](tutorial_det_dl.svg)](tutorial_det_dl.pdf)

# > In the configuration file shown above, the position of the PN junction boundary `pn_r` was calculated as the position where the density is zero.

# ## Display the impurity profile defined

r_list = 0u"mm":0.01u"mm":det_r
imp_list = map(r -> let pt::CylindricalPoint{T} = CylindricalPoint(r, 0u"°", z_draw)
    SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, pt) * 1e-6u"cm^-3"
end, r_list)
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "r / mm", ylabel = "Impurity density / cm⁻³",
          limits = ((0, 10), (-2e10, 1e11)))
lines!(ax, ustrip.(u"mm", r_list), ustrip.(u"cm^-3", imp_list), color = :darkblue, linewidth = 2)
vlines!(ax, [ustrip(u"mm", pn_r)], linewidth = 2, linestyle = :dash, color = :darkred,
        label = "PN junction boundary")
axislegend(ax)
fig
#jl save("tutorial_imp_dl.pdf", current_figure()) # hide
#md save("tutorial_imp_dl.pdf", current_figure()) # hide
#md save("tutorial_imp_dl.svg", current_figure()); nothing # hide
#md # [![tutorial_imp_dl](tutorial_imp_dl.svg)](tutorial_imp_dl.pdf)

# ## Display the mobility curve
using SolidStateDetectors: Electron, Hole
cdm = sim.detector.semiconductor.charge_drift_model
depth_list = 0u"mm":0.01u"mm":(det_r-pn_r)
mobility_list = map(depth -> let pt::CartesianPoint{T} = CartesianPoint(det_r - depth, 0, z_draw)
    (µe = SolidStateDetectors.calculate_mobility(cdm, pt, Electron) * 10000u"cm^2/(V*s)",
     µh = SolidStateDetectors.calculate_mobility(cdm, pt, Hole) * 10000u"cm^2/(V*s)")
end, depth_list)
depths = ustrip.(u"mm", depth_list)
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Depth to surface / mm", ylabel = "Mobility / cm² / Vs",
          limits = ((0, ustrip(u"mm", det_r - pn_r)), nothing))
lines!(ax, depths, ustrip.(getfield.(mobility_list, :µh)), label = "Hole", linewidth = 4)
lines!(ax, depths, ustrip.(getfield.(mobility_list, :µe)), label = "Electron", linewidth = 4)
axislegend(ax, position = :lt)
fig
#jl save("tutorial_mob_dl.pdf", current_figure()) # hide
#md save("tutorial_mob_dl.pdf", current_figure()) # hide
#md save("tutorial_mob_dl.svg", current_figure()); nothing # hide
#md # [![tutorial_mob_dl](tutorial_mob_dl.svg)](tutorial_mob_dl.pdf)

# ## Calculate the electric field
# To simulate the inactive layer, we need to calculate the electric field on a very fine grid.
calculate_electric_potential!(sim, max_n_iterations = 10, grid = Grid(sim), verbose = false, depletion_handling = true)
g = sim.electric_potential.grid
ax1, ax2, ax3 = g.axes
bulk_tick_dis = 0.05u"mm"
dl_tick_dis   = 0.01u"mm"
user_additional_ticks_ax1 = sort(vcat(ax1.interval.left*u"m":bulk_tick_dis:pn_r, pn_r:dl_tick_dis:ax1.interval.right*u"m"))
user_ax1 = typeof(ax1)(ax1.interval, SolidStateDetectors.to_internal_units.(user_additional_ticks_ax1))
user_g = typeof(g)((user_ax1, ax2, ax3))
calculate_electric_potential!(sim, refinement_limits = 0.1, grid = user_g, depletion_handling = true)
calculate_electric_field!(sim)
calculate_weighting_potential!(sim, 1, depletion_handling = true)
calculate_weighting_potential!(sim, 2, depletion_handling = true);
# The same dashed marker goes on all four panels, so it is a function of the axis.
mark_pn!(ax) = vlines!(ax, [ustrip(u"mm", pn_r)], linewidth = 2, linestyle = :dash,
                       color = :darkred, label = "PN junction boundary")

fig = Figure(size = (900, 700))
ap = scalarpotentialplot(fig[1, 1], sim.imp_scale; φ = 0, axis = (; title = "impurity scale"))
mark_pn!(ap.axis)
ap = scalarpotentialplot(fig[1, 2], sim.point_types; φ = 0, axis = (; title = "point types"))
mark_pn!(ap.axis)
ap = scalarpotentialplot(fig[2, 1], sim.electric_potential; axis = (; title = "electric potential"))
mark_pn!(ap.axis)
ap = electricfieldplot(fig[2, 2], sim.electric_field; colorrange = (0, 100 * 2000),
                       axis = (; title = "electric field"))
mark_pn!(ap.axis)
fig
#jl save("tutorial_ef_dl.pdf", current_figure()) # hide
#md save("tutorial_ef_dl.pdf", current_figure()) # hide
#md save("tutorial_ef_dl.svg", current_figure()); nothing # hide
#md # [![tutorial_ef_dl](tutorial_ef_dl.svg)](tutorial_ef_dl.pdf)

# ## Pulse shape and CCE simulation
depth_list = 0.1u"mm":0.1u"mm":(det_r-pn_r)
totTime = 5u"µs"
totEnergy = 1u"keV" # --> simulating ~339 carrier pairs
N = Int(totEnergy÷2.95u"eV")

pulse_fig = Figure(size = (1000, 400))
pulse_ax = Axis(pulse_fig[1, 1], xlabel = "Time / ns", ylabel = "Amplitude / e")
eff_list = map(depth -> begin
    r = det_r-depth
    energy_depos = fill(2.95u"eV", N)
    starting_positions = repeat([CartesianPoint(r, 0, z_draw)], N)
    evt = Event(starting_positions, energy_depos)
    simulate!(evt, sim, Δt = 1u"ns", max_nsteps = round(Int, totTime / 1u"ns"), diffusion = true, end_drift_when_no_field = false, self_repulsion = false)
    pulse = evt.waveforms[1]
    lines!(pulse_ax, ustrip.(u"ns", pulse.time), ustrip.(pulse.signal), linewidth = 2,
           label = "Depth: $(round(typeof(depth), depth, digits = 1))")
    maximum(pulse.signal)/N
end, depth_list)
axislegend(pulse_ax, position = :rt)

cce_ax = Axis(pulse_fig[1, 2], xlabel = "Depth to surface / mm",
              ylabel = "Charge collection efficiency")
# `eff_list` inherits the signal's `e` from `maximum(pulse.signal)/N`, and a ratio labelled "1.0e"
# reads as a charge. This is what `unitformat = :nounit` did on the Plots side.
lines!(cce_ax, ustrip.(u"mm", depth_list), ustrip.(eff_list), linewidth = 2, color = :black)
pulse_fig
#jl save("tutorial_pulse_cce_dl.pdf", current_figure()) # hide
#md save("tutorial_pulse_cce_dl.pdf", current_figure()) # hide
#md save("tutorial_pulse_cce_dl.svg", current_figure()); nothing # hide
#md # [![tutorial_pulse_cce_dl](tutorial_pulse_cce_dl.svg)](tutorial_pulse_cce_dl.pdf)
