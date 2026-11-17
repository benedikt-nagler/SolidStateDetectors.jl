# # Simulation Chain: Inverted Coax Detector

using CairoMakie
using SolidStateDetectors
using Unitful

T = Float32
sim = Simulation{T}(SSD_examples[:InvertedCoax])

# `:wireframe` rather than the default `:csg`: CairoMakie draws a translucent mesh triangle by
# triangle, so a 3-D detector reads as a spiky web (see the warning in the Plotting manual).
detectorplot(sim.detector; seriestype = :wireframe, figure = (; size = (700, 700)))
#jl save("tutorial_det.pdf", current_figure()) # hide
#md save("tutorial_det.pdf", current_figure()) # hide
#md save("tutorial_det.svg", current_figure()); nothing # hide
#md # [![tutorial_det](tutorial_det.svg)](tutorial_det.pdf)

# By default only the contacts and passives are shown in detector plots. One can also take a look at the semiconductor. The following plots show the semiconductor and slices in the y and z direction. 

# One object in three representations, so the legend would say "Semiconductor" three times and eat
# a third of the width; `legend = false` gives the space back to the plots.
fig = Figure(size = (1000, 350))
detectorplot(fig[1, 1], sim.detector.semiconductor; seriestype = :samplesurface, legend = false)
detectorplot(fig[1, 2], sim.detector.semiconductor; seriestype = :slice, y = 0u"mm",
             markersize = 2, legend = false)
detectorplot(fig[1, 3], sim.detector.semiconductor; seriestype = :slice, z = 40u"mm",
             markersize = 2, legend = false)
fig
#jl save("tutorial_semiconductor.pdf", current_figure()) # hide
#md save("tutorial_semiconductor.pdf", current_figure()) # hide
#md save("tutorial_semiconductor.svg", current_figure()); nothing # hide
#md # [![tutorial_semiconductor](tutorial_semiconductor.svg)](tutorial_semiconductor.pdf)

# One can also have a look at how the initial conditions look like on the grid (it starts with a very coarse grid):

apply_initial_state!(sim, ElectricPotential) # optional
fig = Figure(size = (1000, 800))
scalarpotentialplot(fig[1, 1], sim.electric_potential) # initial electric potential (boundary conditions)
scalarpotentialplot(fig[1, 2], sim.point_types) # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
scalarpotentialplot(fig[2, 1], sim.q_eff_imp) # charge density distribution
scalarpotentialplot(fig[2, 2], sim.ϵ_r) # dielectric distribution
fig
#jl save("tutorial_initial_condition.pdf", current_figure()) # hide
#md save("tutorial_initial_condition.pdf", current_figure()) # hide
#md save("tutorial_initial_condition.svg", current_figure()); nothing # hide
#md # [![tutorial_initial_condition](tutorial_initial_condition.svg)](tutorial_initial_condition.pdf)


# Next, calculate the electric potential:

calculate_electric_potential!( sim,
                               refinement_limits = [0.2, 0.1, 0.05, 0.01])

fig = Figure(size = (1000, 800))
scalarpotentialplot(fig[1, 1], sim.electric_potential; φ = 20) # calculated electric potential
scalarpotentialplot(fig[1, 2], sim.point_types) # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
scalarpotentialplot(fig[2, 1], sim.q_eff_imp) # charge density distribution
scalarpotentialplot(fig[2, 2], sim.ϵ_r) # dielectric distribution
fig
#jl save("tutorial_calculated_potential.pdf", current_figure()) # hide
#md save("tutorial_calculated_potential.pdf", current_figure()) # hide
#md save("tutorial_calculated_potential.svg", current_figure()); nothing # hide
#md # [![tutorial_calculated_potential](tutorial_calculated_potential.svg)](tutorial_calculated_potential.pdf)

# SolidStateDetectors.jl supports active (i.e. depleted) volume calculation:

is_depleted(sim.point_types)

# 

get_active_volume(sim.point_types) # approximation (sum of the volume of cells marked as depleted)



# ## Partially depleted detectors

# SolidStateDetectors.jl can also calculate the electric potential of a partially depleted detector:

sim_undep = deepcopy(sim)
sim_undep.detector = SolidStateDetector(sim_undep.detector, contact_id = 2, contact_potential = 500); # V  <-- Bias Voltage of Mantle

calculate_electric_potential!( sim_undep,
                               depletion_handling = true,
                               convergence_limit = 1e-6,
                               refinement_limits = [0.2, 0.1, 0.05, 0.01],
                               verbose = false)


fig = Figure(size = (1000, 400))
scalarpotentialplot(fig[1, 1], sim_undep.electric_potential)
scalarpotentialplot(fig[1, 2], sim_undep.point_types)
scalarpotentialplot(fig[1, 3], sim_undep.imp_scale)
fig
#jl save("tutorial_calculated_potential_undep.pdf", current_figure()) # hide
#md save("tutorial_calculated_potential_undep.pdf", current_figure()) # hide
#md save("tutorial_calculated_potential_undep.svg", current_figure()); nothing # hide
#md # [![tutorial_calculated_potential_undep](tutorial_calculated_potential_undep.svg)](tutorial_calculated_potential_undep.pdf)

# 

is_depleted(sim_undep.point_types)

# Compare both volumes:

println("Depleted:   ", get_active_volume(sim.point_types))
println("Undepleted: ", get_active_volume(sim_undep.point_types));

# The depletion voltage can also be estimated from the simulation.

estimate_depletion_voltage(sim)


# ## Electric field calculation

# Calculate the electric field of the fully depleted detector, given the already calculated electric potential:


calculate_electric_field!(sim, n_points_in_φ = 72)

fap = electricfieldplot(sim.electric_field; full_det = true, φ = 0.0,
                        figure = (; size = (700, 700)), colorrange = (0, 800_000))
fieldlinesplot!(fap.axis, sim; full_det = true, φ = 0.0)
fap.figure
#jl save("tutorial_electric_field.pdf", current_figure()) # hide
#md save("tutorial_electric_field.pdf", current_figure()) # hide
#md save("tutorial_electric_field.svg", current_figure()); nothing # hide
#md # [![tutorial_electric_field](tutorial_electric_field.svg)](tutorial_electric_field.pdf)


# ## Simulation of charge drifts

# Any charge drift model can be used for the calculation of the electric field. If no model is explicitly given, the `ElectricFieldChargeDriftModel` is used. Other configurations are saved in their configuration files and can be found under:

# `<package_directory>/examples/example_config_files/ADLChargeDriftModel/<config_filename>.yaml.`

# Set the charge drift model of the simulation:

charge_drift_model = ADLChargeDriftModel()
sim.detector = SolidStateDetector(sim.detector, charge_drift_model)

# Now, let's create a "random" multi-site event:

starting_positions = [ CylindricalPoint{T}( 0.020, deg2rad(10), 0.015 ),
                       CylindricalPoint{T}( 0.015, deg2rad(20), 0.045 ),
                       CylindricalPoint{T}( 0.022, deg2rad(35), 0.025 ) ]
energy_depos = T[1460, 609, 1000] * u"keV" # are needed later in the signal generation

evt = Event(starting_positions, energy_depos);

time_step = 5u"ns"
drift_charges!(evt, sim, Δt = time_step)

fap = detectorplot(sim.detector; seriestype = :wireframe, figure = (; size = (700, 700)))
driftpathplot!(fap.axis, evt.drift_paths)
fap.figure
#jl save("tutorial_drift_paths.pdf", current_figure()) # hide
#md save("tutorial_drift_paths.pdf", current_figure()) # hide
#md save("tutorial_drift_paths.svg", current_figure()); nothing # hide
#md # [![tutorial_drift_paths](tutorial_drift_paths.svg)](tutorial_drift_paths.pdf)


# ## Weighting potential calculation

# We need weighting potentials to simulate the detector charge signal induced by drifting charges. We'll calculate the weighting potential for the point contact and the outer shell of the detector:

for contact in sim.detector.contacts
    calculate_weighting_potential!(sim, contact.id, refinement_limits = [0.2, 0.1, 0.05, 0.01], n_points_in_φ = 2, verbose = false)
end

fig = Figure(size = (900, 700))
for (i, wp) in enumerate(sim.weighting_potentials)
    ap = scalarpotentialplot(fig[1, i], wp)
    detectorplot!(ap.axis, sim.detector; seriestype = :slice, φ = 0)
end
fig
#jl save("tutorial_weighting_potentials.pdf", current_figure()) # hide
#md save("tutorial_weighting_potentials.pdf", current_figure()) # hide
#md save("tutorial_weighting_potentials.svg", current_figure()); nothing # hide
#md # [![tutorial_weighting_potentials](tutorial_weighting_potentials.svg)](tutorial_weighting_potentials.pdf)

# ## Detector Capacitance Matrix

# After the weighting potentials are calculated, the detector capacitance matrix can be calculated in the **Maxwell Capacitance Matrix Notation**:

calculate_capacitance_matrix(sim)

# See [Capacitances](@ref) for more information.

# ## Detector waveform generation

# Given an interaction at an arbitrary point in the detector, we can now simulate the charge drift and the time evolution of the charges induced on the contacts (e.g. on the point contact):

simulate!(evt, sim) # drift_charges + signal generation of all channels

waveformplot(evt.waveforms[1]; linewidth = 1.5,
             axis = (; limits = ((0, 1100), nothing), ylabel = "Charge / e",
                       xticklabelsize = 14, yticklabelsize = 14, xlabelsize = 16, ylabelsize = 16))
#jl save("tutorial_waveforms.pdf", current_figure()) # hide
#md save("tutorial_waveforms.pdf", current_figure()) # hide
#md save("tutorial_waveforms.svg", current_figure()); nothing # hide
#md # [![tutorial_waveforms](tutorial_waveforms.svg)](tutorial_waveforms.pdf)

# SolidStateDetectors.jl also allows separating the observed charge signal into the charge induced by the electrons and the charge induced by the holes.

contact_id = 1
simulate!(evt, sim, signal_unit = u"keV")
ehcontributionplot(evt, sim, contact_id;
                   axis = (; limits = ((0, 1100), nothing),
                             xticklabelsize = 14, yticklabelsize = 14,
                             xlabelsize = 16, ylabelsize = 16))
#jl save("tutorial_waveform_contributions.pdf", current_figure()) # hide
#md save("tutorial_waveform_contributions.pdf", current_figure()) # hide
#md save("tutorial_waveform_contributions.svg", current_figure()); nothing # hide
#md # [![tutorial_waveform_contributions](tutorial_waveform_contributions.svg)](tutorial_waveform_contributions.pdf)