# Plotting

In order to plot geometries or simulation results, load a [Makie.jl](https://docs.makie.org) backend.
[CairoMakie.jl](https://docs.makie.org/stable/explanations/backends/cairomakie) writes publication-quality
vector graphics and needs no display; [GLMakie.jl](https://docs.makie.org/stable/explanations/backends/glmakie)
opens an interactive window and renders three-dimensional geometry better.

First, we have to calculate something in order to demonstrate the plotting tools:
````@setup tutorial
using CairoMakie
using SolidStateDetectors
using Unitful

T = Float32
sim = Simulation{T}(SSD_examples[:InvertedCoax]);
simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01]);
````
````julia
using CairoMakie
using SolidStateDetectors
using Unitful

T = Float32
sim = Simulation{T}(SSD_examples[:InvertedCoax]);
simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01]);
````

Every object below has a named plotting function, and `plot(object)` routes to the same place. The
named form is used here because that is where the keywords are documented.

Each function comes in three forms:

* `detectorplot(det)` builds its own `Figure`, with the title, axis labels and legend set;
* `detectorplot(fig[1, 2], det)` draws into one cell of a figure you already have, which is how the
  multi-panel figures below are made;
* `detectorplot!(ax, det)` draws into an existing axis and sets nothing, which is how one plot is
  laid over another.

Figure-wide settings go through `figure`, axis settings through `axis`: `figure = (; size = (500, 500))`
replaces Plots' `size = (500, 500)`, and `axis = (; title = "…")` sets a title.

Coordinates are drawn in millimetres and angles in degrees throughout, with the unit in the axis
label. Under Plots the units had to be primed with a bare `plot(u"cm", u"cm", u"cm")` call before
plotting; that is no longer needed.

## Detector Plots

The geometry of a detector, together with its environment, can be simply plotted via
````@example tutorial
det = sim.detector
detectorplot(det, figure = (; size = (500, 500)))
````
### Plot Styles

The style of the detector plot can be controlled via the `seriestype` keyword. There are 4 available styles for 3D plotting:

````@example tutorial
fig = Figure(size = (800, 220))
for (i, st) in enumerate((:csg, :wireframe, :mesh3d, :samplesurface))
      extra = st === :samplesurface ? (; markersize = 1, n_samples = 50) : (;)
      ap = detectorplot(fig[1, i], det; seriestype = st, legend = false,
                        axis = (; title = ":$st"), extra...)
      hidedecorations!(ap.axis)
end
fig
````

To generate a 3D plot the `seriestype` can be set to `:csg` (default), `:wireframe`, `:mesh3d`, or `:samplesurface`. `:csg` plots a wireframe on top of a mesh (with no mesh grid lines). For fastest plotting use either `:wireframe` or `:mesh3d` and consider changing `n_arc` (see [Optional Keywords](@ref)). 

!!! warning "Translucent meshes under CairoMakie"
      `:csg` and `:mesh3d` draw the detector as a translucent mesh, and CairoMakie renders one by
      drawing every triangle: the body reads as a dense web of lines rather than a smooth solid.
      There is no attribute that fixes it. For a publication figure use GLMakie, which composites
      transparency correctly, or one of the two styles that do not rely on it, `:wireframe` and
      `:samplesurface`.

!!! tip
      `:csg`, `:wireframe`, and `:mesh3d` are all mesh-based. For geometries containing differences or intersections the recommended `seriestype` is `:samplesurface` which is marker-based. The marker density is set low by default for speed. For increased plot fidelity use `n_samples` above `100` and adjust the `markersize`.

Plotting the detector's semiconductor helps illustrate the difference between each `seriestype`.
````@example tutorial
fig = Figure(size = (800, 220))
for (i, st) in enumerate((:csg, :wireframe, :mesh3d, :samplesurface))
      extra = st === :samplesurface ? (; markersize = 1, n_samples = 50) : (;)
      ap = detectorplot(fig[1, i], det.semiconductor; seriestype = st, legend = false,
                        axis = (; title = ":$st"), extra...)
      hidedecorations!(ap.axis)
end
fig
````

!!! note
    When using mesh-based `seriestype`'s (`:csg`, `:wireframe`, `:mesh3d`), plots are produced by plotting whole primitives. Thus, usually, nicer plots are produced if the geometry consists only of unions of primitives and not differences or intersections. If the geometry contains differences, the resulting negative geometries are plotted with thinner wireframe lines and/or with semi-transparent white mesh faces depending on the `seriestype` used. 

For 2D plotting, use `seriestype = :slice`. Like `:samplesurface`, `:slice` is marker based and will respond to all marker keywords like `markersize`. As the `seriestype` implies, a slice of the detector is plotted. To select the slice value use the `x`, `y` or `z` or `φ` keywords.

````@example tutorial
fig = Figure(size = (900, 240))
detectorplot(fig[1, 1], det; seriestype = :slice, x = 2u"mm", legend = false)
detectorplot(fig[1, 2], det; seriestype = :slice, y = 20u"mm", color = :black, legend = false)
detectorplot(fig[1, 3], det; seriestype = :slice, z = 0u"mm", legend = false)
detectorplot(fig[1, 4], det; seriestype = :slice, φ = 90u"°", markersize = 3, legend = false)
fig
````  
`seriestype = :slice` has all the advantages of `:samplesurface` but with much greater fidelity and speed.
!!! tip
      For `seriestype = :slice`, `markersize` and `linewidth` can be used interchangeably. For very small `markersize`, consider increasing `n_samples`. This `seriestype` is particularly useful for showing contacts over potentials. See examples in [Scalar Potential Plots](@ref)

### How does the plot recipe work?

The detector consists of a semiconductor, `det.semiconductor`, its contacts, `det.contacts`,
and its surrounding objects, `det.passives`.
By default (see [Optional Keywords](@ref)), only the contacts and the surrounding objects are plotted
and the color of the contact primitives is defined internally through their `id`.

Components can also be plotted individually for enhanced style handling. Make the axis once, then
draw each component into it with the mutating form.
````@example tutorial
fig = Figure(size = (500, 500))
ax = Axis3(fig[1, 1], xlabel = "x / mm", ylabel = "y / mm", zlabel = "z / mm",
           aspect = :data, azimuth = 0.7π, elevation = 0.2π)
csgplot!(ax, det.semiconductor.geometry; seriestype = :samplesurface, n_samples = 100, markersize = 2)
csgplot!(ax, det.contacts[1].geometry; seriestype = :mesh3d, linewidth = 0.5, color = :white)
csgplot!(ax, det.contacts[2].geometry; seriestype = :wireframe, n_vert_lines = 5)
fig
````
The camera of an `Axis3` is set with `azimuth` and `elevation`, in radians, where Plots took a
`camera = (40, 55)` pair in degrees.

### Optional Keywords
* `show_semiconductor`: Will display the semiconductor if set to `true`. Default is `false`.
* `show_passives`: Will display the objects surrounding the detector if set to `true`. Default is `true`.
* `seriestype`: Can be `:csg` (default), `:wireframe`, `:mesh3d`, `:samplesurface` or `:slice`.
* `linewidth`: Sets the line width of the edges of the mesh gridlines when using `seriestype = :mesh3d`. When using `seriestype = :csg` or `seriestype = :wireframe`, `linewidth` sets the line width of the wireframe. For `seriestype = :slice`, `linewidth` sets the `markersize`.
* `linecolor`: Sets the line color of the edges of the mesh gridlines when using `seriestype = :mesh3d`. When using `seriestype = :csg` or `seriestype = :wireframe`, `linecolor` sets the line color of the wireframe.
* `color`: Sets the face color of all faces of the mesh, and the marker color of the marker-based styles.
* `alpha`: Sets the opacity of the faces and markers.
* `markersize`: Sets the marker size. The default `markersize` is `3` and `1` for `:samplesurface` and `:slice` respectively.
* `legend`: Set to `false` to leave the legend off.
* `n_arc`: Controls the discretization of curved objects in a mesh. Each full ellipse is divided into `n_arc` segments. Partial ellipses are drawn with a proportional number (`n_arc*f` with `f<1`) of segments.`n_arc = 40` is the default value. Smaller `n_arc` values will result in faster plotting, specially if the geometry contains tori or ellipsoids.
* `n_vert_lines`: Controls the number of wireframe "vertical" lines in a mesh. `n_vert_lines = 2` is the default value. A maximum of `n_arc*f` vertical lines can be drawn on each curved object. This keyword is ignored by polygons.
* `n_samples`: Controls the marker density. `n_samples = 40` is the default value for `:samplesurface`. For `:slice` the default value is set dynamically to `100` or `200`. Reduce this value for faster plotting. Consider increasing `markersize` and/or `markeralpha` to compensate for the visual impact of lower marker densities. Note that the marker density is intended to be even across all dimensions. Therefore, visual distortions will occur if the aspect ratio of the axes is far from unity.
* `x,y,z,φ`: In combination with `seriestype = :slice` selects the slice to be plotted.
* `full_det`: For `φ` slices determines if the cross section in `φ` is extended to "negative" `r`. Default value is `false`.
* `projection`: Determines if `z` slices will be shown in a `:polar` projection. Default value is `:none`.


## Scalar Potential Plots

When calculating the electric potential, a set of scalar potentials are stored in the `Simulation` object, e.g. `sim.electric_potential`, `sim.point_types`, and `sim.q_eff_imp`. SolidStateDetectors.jl has plot recipes for these scalar potentials.

````@example tutorial
fap = scalarpotentialplot(sim.electric_potential; φ = 0)
detectorplot!(fap.axis, sim.detector; seriestype = :slice, φ = 0)
fap.figure
````

There are several keyword arguments that can be passed to adjust the potential plots.
First, the simulated potentials are simulated on a three-dimensional `Grid`, but the plots will always be a two-dimensional cross section (for a `CylindricalGrid` at constant `r`, `φ` or `z`, for a `CartesianGrid3D` at constant `x`, `y` or `z`).
The cross section can be specified in the `plot` command. A cross section of the `ElectricPotential` at `z = 30mm` can be plotted via
````@example tutorial
fap = scalarpotentialplot(sim.electric_potential; z = 30u"mm")
detectorplot!(fap.axis, det; seriestype = :slice, z = 30u"mm", projection = :polar, markersize = 4)
fap.figure
````
!!! warning
      Only plot the detector over scalar potentials with `seriestype = :slice`. When plotting the detector over a polar projection of the potential (such as above), `proj = :polar` must be selected in the detector plot.

In addition to the usual Makie attributes, plots of scalar potentials can take two additional keyword arguments, `contours_equal_potential` and `full_det`, which both are of type `Bool`:
* `contours_equal_potential`: If this is set to `true`, the plot will additionally display contour lines at points with equal potential value.
* `full_det`: By default, cross sections in `φ` are only displayed for positive radii `r`. If `full_det` is set to `true`, cross sections in `φ` are also extended to "negative" `r`, by additionally plotting the potential values at `φ + 180°` on the left side of the plot.

````@example tutorial
scalarpotentialplot(sim.electric_potential; φ = 30u"°",
                    contours_equal_potential = true, full_det = true,
                    contour_levels = 34)
````


## Electric Field Plots

The magnitude of the electric field can be plotted using the same syntax as [Scalar Potential Plots](@ref).

````@example tutorial
electricfieldplot(sim.electric_field; full_det = true, colorrange = (0, 5e5))
````

In addition, SolidStateDetectors.jl offers the possibilities to plot electric field lines via [`fieldlinesplot`](@ref).
This is done by spawning charges close to the surface of the contacts and simulating their drift
parallel to the electric field, using [`ElectricFieldChargeDriftModel`](@ref).

The lines are drawn in white, to be laid over a field-magnitude plot:
````@example tutorial
fap = electricfieldplot(sim.electric_field; full_det = true, colorrange = (0, 5e5))
fieldlinesplot!(fap.axis, sim; full_det = true)
fap.figure
````

In addition to the usual Makie attributes and the syntax for specifying the cross section and
`full_det` (see [Scalar Potential Plots](@ref)), the field-line plot understands the following
arguments.
* `sampling`: Specifies the steps at which the contacts are sampled to generate equally spaced charges at the surface. The default is `2u"mm"`, but the optimal value depends on the geometry of the detector and contacts. If no unit is given, `sampling` is interpreted in units of meter.
* `offset`: The charges are created on the surface and have to be moved slightly inside the semiconductor to be able to drift (charges that are in the contacts will not drift). This keyword defines how much the charges will be moved inside along the normal vector of the surface. The default is `0.5u"mm"`, but the optimal value again depends on the detector geometry. Unitless quantities are interpreted in units of meter.
* `skip_contact`: Detectors will usually have positively and negatively biased contacts. The charges should only be spawned on either of them, but not both. This keyword defines which contact should be skipped in the charge spawn algorithm. Default is `1`, i.e. that no charges will be spawned at the surface of contact with `id = 1`.
* `max_nsteps`: The generation of electric field lines is based on the charge drift code, which requires a maximum number of steps, which is set to `5000` by default here. If the drift ends before the charges reach a contact, the electric field line will not be fully displayed. In this case, it is recommended to increase `max_nsteps`.
* `Δt`: The charges are drifted with a normalized mobility of $1\,\text{cm}^2\text{/Vs}$ (see [`ElectricFieldChargeDriftModel`](@ref)). If the individual steps are too large, it is recommended to reduce `Δt`.



In case the electric field line plots do not look good, adjust `sampling`, `offset` and `max_nsteps` until obtaining the desired result.
````@example tutorial
fap = electricfieldplot(sim.electric_field; full_det = true, colorrange = (0, 5e5),
                        figure = (; size = (500, 500)))
fieldlinesplot!(fap.axis, sim; full_det = true, sampling = 3u"mm", offset = 2u"mm")
fap.figure
````

## Event Plots

SolidStateDetectors.jl also provides plot recipes to display charge drift paths resulting from simultaneous energy deposits
in the semiconductor body of the detector. To simulate and plot an individual event, the `Event` struct is recommended.

````@example tutorial
evt = Event([CartesianPoint{T}(0.01,0.01,0.075)], [2u"MeV"])
simulate!(evt, sim)
fap = detectorplot(sim.detector, figure = (; size = (500, 500)))
driftpathplot!(fap.axis, evt.drift_paths)
fap.figure
````

The trajectories take the usual line attributes:
````@example tutorial
fap = detectorplot(sim.detector, figure = (; size = (500, 500)))
driftpathplot!(fap.axis, evt.drift_paths; linewidth = 2, linestyle = :dash, markersize = 6)
fap.figure
````

There are also plot recipes for plotting the simulated waveforms of the `Event`:
````@example tutorial
waveformplot(evt.waveforms; labels = ["Contact 1", "Contact 2"])
````

The length of the waveforms is given by the length of the charge drift.
By default, no baseline and no tail are added to the waveforms. However, this might be desired in waveforms plots.
The waveforms can be extended by calling [`add_baseline_and_extend_tail`](@ref) on the waveforms.
````@example tutorial
waveformplot(add_baseline_and_extend_tail.(evt.waveforms, 0, 400);
             linewidth = 4, linestyle = :dash)
````
