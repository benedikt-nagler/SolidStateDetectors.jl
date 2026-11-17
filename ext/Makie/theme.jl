# The SolidStateDetectors colour scheme.
#
# The Plots recipes colour a contact by `seriescolor --> contact.id`, i.e. "entry `id` of the current
# palette". Makie's equivalent, `Cycled(i)`, indexes the *Wong* palette, which is a different set of
# colours — so a naive port silently recolours every detector. Worse, both are *defaults*: they move
# whenever either package changes its palette.
#
# So the palette is pinned here instead. The values are Plots' default palette, read out of
# `Plots.palette(:default)` and cross-checked against the committed golden specs (contact 1 is
# `#009afa`, contact 2 `#e36f47` in `test/parity/golden/detector_cyl_*.json`). Contacts therefore
# keep exactly the colours they have today, and keep them regardless of what either backend does to
# its defaults later.

"""
    SSD_PALETTE

Colours indexed by contact `id`, pinned to Plots' default palette so contact colours are identical
across the two backends and stable across backend versions.
"""
const SSD_PALETTE = Makie.RGBAf[
    Makie.to_color("#009afa"), Makie.to_color("#e36f47"), Makie.to_color("#3ea44e"),
    Makie.to_color("#c371d2"), Makie.to_color("#ac8e18"), Makie.to_color("#00aaae"),
    Makie.to_color("#ed5e93"), Makie.to_color("#c68225"), Makie.to_color("#00a98d"),
    Makie.to_color("#8e971e"),
]

"""
    contact_color(id) -> RGBAf

The colour of the contact with this `id`, wrapping around for detectors with more contacts than the
palette has entries (which is what Plots' palette cycling does too).
"""
contact_color(id::Integer) = SSD_PALETTE[mod1(Int(id), length(SSD_PALETTE))]

"""
    color_or_default(given, fallback)

`given` unless it is `Makie.automatic`, in which case `fallback`.

Every recipe here takes its colours as `automatic` so that the figure-level wrapper can read the
*same* value back when it builds the legend: a legend that describes one colour while the plot draws
another is the failure mode explicit legend elements otherwise invite (plan §1.3).
"""
color_or_default(given, fallback) = given === Makie.automatic ? fallback : given

"""
    SSD_COLORS

Fixed colours for the parts that are not palette-indexed. Same values the Plots recipes use
(`:grey`, `:silver`, `:red`, `:green`).
"""
const SSD_COLORS = (
    semiconductor   = Makie.to_color("#808080"),   # :grey
    passive         = Makie.to_color("#c0c0c0"),   # :silver
    electron        = Makie.to_color("#ff0000"),   # :red
    hole            = Makie.to_color("#008000"),   # :green
    fieldline       = Makie.to_color("#ffffff"),   # :white, over a field-magnitude plot
    fieldline_alone = Makie.to_color("#1a1a1a"),   # on an empty axis
    charge_cloud    = Makie.to_color("#0000ff"),   # :blue
    particle_source = Makie.to_color("#808080"),   # :gray, outlined in black
    emission        = Makie.to_color("#008000"),   # :green, the ray or cone a source emits into
)

"""
    SSD_ALPHAS

Opacities of the geometry layers, matching the `fillalpha` values of the Plots recipes: contacts and
passives are drawn translucent so the semiconductor behind them stays visible, while a legend swatch
is always drawn opaque (which is what the empty-`:path`-series hack in the Plots recipes existed to
achieve — see plan §1.3).

`semiconductor` is `0.5`, not `0.2` like the other two: the Plots `Semiconductor` recipe is the one
that sets no `fillalpha` at all, so it falls through to the `Mesh` recipe's own default. Its golden
spec records the difference (`detector_cyl_semiconductor`), and it is why a semiconductor reads as
the solid body with the contacts as faint shells on it.
"""
const SSD_ALPHAS = (
    contact       = 0.2,
    passive       = 0.2,
    semiconductor = 0.5,
    mesh          = 0.5,
    subtracted    = 0.2,
    samplesurface = 0.2,
    legend        = 1.0,
)

"""
    SSD_LINEWIDTHS

Line widths of the geometry representations. Makie measures these in px at 72 dpi and Plots does not
(plan §2), so these are the one place the two scales are reconciled; they are calibrated against the
gallery rather than copied from the Plots recipes.

`wireframe` is the black outline of `:csg` and `:wireframe`, `meshgrid` the lighter cell grid
`:mesh3d` draws over its faces, and `subtracted` the thinner outline a primitive removed by a
difference is drawn with. `driftpath` and `charge_cloud` are the two the Plots recipes state outright
(3 and 1), kept as they are because a line width in Plots is already close to Makie's px; the two
waveform widths are Plots' 1 and 2 kept the same way.
"""
const SSD_LINEWIDTHS = (
    wireframe    = 1.0,
    meshgrid     = 0.5,
    subtracted   = 0.8,
    driftpath    = 3.0,
    charge_cloud = 1.0,
    waveform     = 1.0,
    contribution = 2.0,
    fieldline    = 1.0,
)

"""
    SSD_MARKERSIZES

Marker sizes of the marker-based representations, in Makie px.

`point_scale` converts the adaptive marker size the Plots `AbstractVector{<:CartesianPoint}` recipe
computes (which is backend-agnostic and kept) from Plots' marker scale to Makie's. The four sizes
below it are the corresponding Plots values already multiplied by it, so a drift-path end point stays
the same fraction of a charge-cloud center that it was; the comment gives the Plots value each came
from.
"""
const SSD_MARKERSIZES = (
    samplesurface       = 3.0,
    slice               = 2.5,
    point_scale         = 2.8,
    drift_endpoint      = 11.0,   # Plots' 4
    particle_source     = 14.0,   # Plots' 5
    charge_cloud        = 17.0,   # Plots' 6
    charge_cloud_center = 28.0,   # Plots' 10
)

"""
    SSD_COLORMAPS

Colormaps not already carried by `_get_potential_plot_information` (which returns `:viridis` or
`:inferno` per scalar field). `contours` is the `:thermal` gradient the Plots recipes force on the
equipotential contour overlay.
"""
const SSD_COLORMAPS = (
    field    = :inferno,
    contours = :thermal,
)

"""
    ssd_theme(; kwargs...) -> Makie.Theme

The package theme: it pins the palette (D3) and nothing else, so a user's own theme keeps control of
fonts, sizes and backgrounds. Extra `kwargs` are merged in.

Use it as `Makie.with_theme(ssd_theme()) do … end`, or `Makie.set_theme!(ssd_theme())` for a session.
"""
ssd_theme(; kwargs...) = Makie.Theme(
    palette = (color = SSD_PALETTE, patchcolor = SSD_PALETTE),
    Lines = (cycle = [:color],),
    Scatter = (cycle = [:color],),
    kwargs...,
)
