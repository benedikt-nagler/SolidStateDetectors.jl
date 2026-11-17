#=
    This is a script to test the whole simulation chain with all standard detector types.
    The script also produces some output plots but only very basic ones, since it
    is only to test the core functionality of the package. (So not detector specific plots)
=#

outputdir = joinpath(ENV["HOME"], "tmp/test_solidstatedetectors.jl/")
mkpath(outputdir)
@info "Test output dir: $outputdir"

@info "Loading packages"
using CairoMakie
using SolidStateDetectors; SSD = SolidStateDetectors
using LinearAlgebra

T = Float32
@info "Testing now for Float32:"

for key in [:InvertedCoax, :BEGe, :Coax, :CGD, :Hexagon]

    @info "Now test detector type: $key"

    sim = Simulation(SSD_examples[key]);
    S = SSD.get_coordinate_system(sim)

    # The three initial fields side by side, in one figure. A Cartesian setup is cut at y = 2 mm;
    # a cylindrical one takes the default cross section.
    cut = S == SSD.Cartesian ? (; y = 0.002) : (;)
    function field_panels(sim)
        fig = Figure(size = (1000, 600))
        for (i, f) in enumerate((sim.electric_potential, sim.q_eff_imp, sim.point_types))
            scalarpotentialplot(fig[1, i], f; cut...)
        end
        fig
    end

    apply_initial_state!(sim, ElectricPotential)
    save(joinpath(outputdir, "$(key)_0_init_setup.png"), field_panels(sim))

    nrefs = if key in [:InvertedCoax, :Heaxgon, :CGD]
        0:3
    else
        0:1
    end
    for nref in nrefs
        SSD.update_till_convergence!(sim, ElectricPotential)
        save(joinpath(outputdir, "$(key)_1_Electric_Potential_$(nref)_refinements.png"),
             field_panels(sim))
        if nref != nrefs[end]
            SSD.refine!(sim, ElectricPotential, (100, 100, 100), (1e-4, 1e-4, 1e-4), update_other_fields = true)
        end
        @show size(sim.electric_potential.grid)
    end

    for contact in sim.detector.contacts
        calculate_weighting_potential!(sim, contact.id, refinement_limits = key == :Coax ? missing : [0.2], verbose = true)
    end
    wpot_fig = Figure(size = (1000, 1000))
    for (i, contact) in enumerate(sim.detector.contacts)
        scalarpotentialplot(wpot_fig[1, i], sim.weighting_potentials[contact.id]; cut...)
    end
    save(joinpath(outputdir, "$(key)_2_Weighting_Potentials.png"), wpot_fig)

    calculate_electric_field!(sim)


    slice = S == SSD.Cylindrical ? (; φ = 30) : (; y = 0)
    fap = electricfieldplot(sim.electric_field; slice...)
    save(joinpath(outputdir, "$(key)_3_Electric_Field_strength.png"), fap.figure)
    fieldlinesplot!(fap.axis, sim; slice...)
    save(joinpath(outputdir, "$(key)_3_1_Electric_Field_Lines.png"), fap.figure)

    sim.detector = SolidStateDetector(sim.detector, ADLChargeDriftModel(T = T))

    pos = if key == :InvertedCoax
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.02, deg2rad(10), 0.025 ) ]
    elseif key == :CGD
        CartesianPoint{T}[ CartesianPoint{T}( 0.003, 0.003, 0.003 ) ]
    elseif key == :BEGe
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.016, deg2rad(10), 0.015  ) ]
    elseif key == :Coax
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.016, deg2rad(10), 0.005  ) ]
    elseif key == :Spherical
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.00, deg2rad(0), 0.0  ) ]
    elseif key == :Hexagon 
        CartesianPoint{T}[ CartesianPoint{T}( 0.0004, 0.0, 0.0 )]
    end
    energy_depos = T[1460]
    @assert in(pos[1], sim.detector) "Test point $(pos[1]) not inside the detector $(key)."

    event = Event(pos, energy_depos)
    simulate!(event, sim)

    drift = detectorplot(sim.detector)
    driftpathplot!(drift.axis, event.drift_paths)
    save(joinpath(outputdir, "$(key)_4_charge_drift.png"), drift.figure)

    save(joinpath(outputdir, "$(key)_5_induced_signals.png"),
         waveformplot([event.waveforms...]; linewidth = 1.5,
                      figure = (; size = (1200, 600))).figure)

end

@info "Finished testing."
@info "Test output saved in: $outputdir"
