# Tier-1 gate: the rules that decide whether a Makie `CaseSpec` matches the Plots golden spec.
#
# The two specs are *not* expected to be identical — the two backends record the same picture
# differently, and a few of the differences are deliberate (plan A1). This file is where every one of
# those bridging rules is written down once, instead of being re-derived in each testset:
#
#   * **Units.** Plots strips units onto the axis (guide `"r"`, unit `"m"`) and leaves the data in
#     internal units; the Makie recipes rescale their data to display units and write the unit into
#     the label (D2), so the expected label is `"r / mm"` and the expected limits are the golden ones
#     × 1000. `display_unit` — pinned by its own unit tests — is what makes the two commensurable.
#   * **Orientation.** A Plots heatmap `z` is `(length(y), length(x))`, a Makie one
#     `(length(x), length(y))`, so the data comparison undoes the transpose (plan §1.7).
#   * **Cells.** Plots was handed cell centers, Makie is handed cell edges, so the `x`/`y` *data* of
#     a heatmap is not comparable at all; the axis limits are, and they are what the picture shows.
#   * **F3.** A Cartesian contour plot loses its axis guides under Plots (plan §3.5). Keeping them is
#     an accepted improvement, so an empty golden guide next to a non-empty Makie label is recorded
#     as a deviation rather than a failure.
#   * **Titles.** Plots printed the cut position in the unit the caller passed (`"@ z = 0.0375 m"`),
#     the Makie wrapper prints it in the display unit (`"@ z = 37.5 mm"`) so it agrees with the axis
#     labels, so the golden title is converted before comparing (`_expected_title`).
#
# This file is plain code, not a module: it is included into a scope that already has
# `PlotParitySpec` (via `MakieSpec.PlotParitySpec`), `Unitful` and the extension module.

# Numbers in a golden spec are stored to 8 significant digits (`PlotParitySpec._json`), so no
# comparison against one may be exact.
const SPEC_ROUNDING_RTOL = 1e-7

"""
    case_diffs(golden, actual; ext) -> (diffs, deviations)

Compare a Makie `CaseSpec` against the Plots golden spec of the same case.

`diffs` is one line per parity failure — empty means the case passes. `deviations` lists the
accepted, deliberate differences that were exercised, so a testset can report them instead of
silently swallowing them. `ext` is the loaded `SolidStateDetectorsMakieExt` module, which supplies
the display-unit convention both sides have to agree on.
"""
function case_diffs(golden::CaseSpec, actual::CaseSpec; ext)
    diffs, deviations = String[], String[]
    isempty(golden.error) || return (["golden spec records an error: $(golden.error)"], deviations)
    isempty(golden.note) || return (["golden spec is not a valid reference: $(golden.note)"], deviations)

    if isempty(actual.axes) || isempty(golden.axes) || length(actual.axes) != length(golden.axes)
        push!(diffs, "axis count $(length(actual.axes)) != $(length(golden.axes))")
        return (diffs, deviations)
    end

    # Two field-line cases legitimately draw *nothing*: the CGD contacts are slabs perpendicular to
    # x, so no contact surface meets an x cross section and no charge is ever spawned. Plots then
    # leaves the subplot with no attributes set at all - no title, no guides, no aspect - so the only
    # comparable property is the emptiness itself, which is the parity property (plan §3.5).
    if isempty(golden.series)
        if isempty(actual.series)
            push!(deviations, "empty: golden drew nothing, only emptiness is compared")
        else
            push!(diffs, "golden drew nothing, makie drew $(length(actual.series)) series")
        end
        return (diffs, deviations)
    end

    # -- what only a figure-level wrapper can set (plan §1.1) ------------------------------------
    # Per panel: everything up to and including the axis box. A `GridLayout` case (grids, step 8) has
    # several; every other case in the matrix has one, and then this loop runs once.
    for i in eachindex(golden.axes)
        append!(diffs, _panel_diffs(golden, actual, i; ext, deviations))
    end

    g, a = golden.axes[1], actual.axes[1]

    # -- legend ----------------------------------------------------------------------------------
    # A Makie recipe produces no legend entries at all (plan §1.3), so this checks the figure-level
    # wrapper built them, in the same order and with the same text.
    g.legend_entries == a.legend_entries ||
        push!(diffs, "legend entries $(a.legend_entries) != $(g.legend_entries)")

    # -- colorbar and colour scale ---------------------------------------------------------------
    # `sp[:colorbar]` is `:best` on every Plots subplot, drawn or not, and a 3-D subplot even gets a
    # `colorbar_title` (the `z` axis unit, `"m"`). Neither says a colorbar exists, so what is asked
    # here is whether the golden plot maps values through a colour scale at all.
    if _maps_colors(golden)
        if a.colorbar == "none"
            push!(diffs, "no colorbar (golden title $(repr(g.colorbar_title)))")
        else
            g.colorbar_title == a.colorbar_title ||
                push!(diffs, "colorbar label $(repr(a.colorbar_title)) != $(repr(g.colorbar_title))")
        end
    end
    if g.clims !== nothing
        if a.clims === nothing
            push!(diffs, "colorrange absent (golden clims $(g.clims))")
        elseif !isapprox(a.clims, g.clims; rtol = 1e-4)
            push!(diffs, "colorrange $(a.clims) != $(g.clims)")
        end
    elseif a.clims !== nothing
        # Plots let the backend take the range from the plotted slice; an explicit Makie colorrange
        # is only wrong if it does not cover that slice. The comparison needs a tolerance rather
        # than a bare `<=`: the golden extrema were rounded to 8 significant digits on write
        # (`_json`), so a colorrange that *is* the data range — which is what the field-magnitude
        # recipe sets — sits a rounding step below its own golden value half the time.
        zex = _field_series(golden) === nothing ? nothing : _field_series(golden).z.extrema
        if zex !== nothing
            slack = SPEC_ROUNDING_RTOL * max(abs(zex[1]), abs(zex[2]))
            (a.clims[1] <= zex[1] + slack && a.clims[2] >= zex[2] - slack) ||
                push!(diffs, "colorrange $(a.clims) does not cover the plotted data $(zex)")
        end
    end
    append!(diffs, compare_colormap(g.colormap, a.colormap))

    # -- the drawn data --------------------------------------------------------------------------
    gh, ah = _field_series(golden), _field_series(actual)
    if gh === nothing && ah === nothing
        # Neither side drew a field, so this is geometry: compare the bounding box of everything
        # drawn, which is the only thing the two backends hold commensurably.
        append!(diffs, _extent_diffs(golden, actual; ext))
    elseif gh === nothing || ah === nothing
        push!(diffs, "field series present on only one side " *
                     "(golden: $(gh !== nothing), makie: $(ah !== nothing))")
    elseif g.projection == "polar"
        append!(diffs, _polar_data_diffs(gh.z, ah.z; ext, deviations))
    else
        append!(diffs, compare_data(gh.z, ah.z; transposed = true, label = "field z"))
    end

    if _has_contours(golden) != _has_contours(actual)
        push!(diffs, "contour overlay: makie $(_has_contours(actual)), golden $(_has_contours(golden))")
    end

    (diffs, deviations)
end

"""
    _panel_diffs(golden, actual, i; ext, deviations) -> Vector{String}

Compare panel `i` of the two specs: title, projection, axis labels, limits and aspect.

Everything here is a property of the `Axis`, which is to say of the figure-level wrapper (plan §1.1).
A case with a `GridLayout` has several panels and this runs once per panel; every other case has one.
"""
function _panel_diffs(golden::CaseSpec, actual::CaseSpec, i::Int; ext, deviations)
    diffs = String[]
    g, a = golden.axes[i], actual.axes[i]
    at = length(golden.axes) == 1 ? "" : "panel $i "

    if _is_grid(golden)
        append!(diffs, _grid_panel_diffs(g, a, i; ext, deviations))
    else
        expected_title = _expected_title(g.title; ext)
        expected_title == a.title ||
            push!(diffs, "$(at)title $(repr(a.title)) != $(repr(expected_title))")
    end

    g.projection == a.projection ||
        push!(diffs, "$(at)projection $(repr(a.projection)) != $(repr(g.projection))")

    if g.projection != "polar"   # a PolarAxis has no labels, no aspect and no Cartesian limits
        _is_grid(golden) || append!(diffs, _label_diffs(g, a, golden, at; ext, deviations))

        for (dim, glims, alims, gunit) in (("x", g.xlims, a.xlims, g.xunit),
                                           ("y", g.ylims, a.ylims, g.yunit))
            glims === nothing && continue
            expected = _expected_lims(glims, gunit; ext)
            if alims === nothing
                push!(diffs, "$(at)$(dim) limits absent (golden $(expected))")
            elseif !isapprox(alims, expected; rtol = 1e-4)
                push!(diffs, "$(at)$(dim) limits $(alims) != $(expected)")
            end
        end

        if _equal_aspect(g.aspect_ratio) && !_equal_aspect(a.aspect_ratio)
            push!(diffs, "$(at)aspect $(repr(a.aspect_ratio)) is not equal-scaled")
        elseif !_equal_aspect(g.aspect_ratio) && _equal_aspect(a.aspect_ratio)
            push!(deviations, "aspect: equal data scaling where Plots left the axis free")
        end
    end
    diffs
end

# The axis labels, with the two accepted losses on the Plots side folded in (F3 and F4).
function _label_diffs(g, a, golden::CaseSpec, at::AbstractString; ext, deviations)
    diffs = String[]
    dims = [("x", g.xguide, g.xunit, a.xguide), ("y", g.yguide, g.yunit, a.yguide)]
    # A third axis exists only in 3-D, and only there does a lost `zguide` mean anything.
    g.projection == "3d" && push!(dims, ("z", g.zguide, g.zunit, a.zguide))
    for (dim, gname, gunit, aname) in dims
        expected = _expected_label(gname, gunit; ext)
        if expected == aname
            # match
        elseif !isempty(gname) && !isempty(aname) && gname == _cut_dimension(g.title)
            push!(deviations, "F5: $(dim)label $(repr(aname)) where Plots named $(repr(gname)), " *
                              "the dimension the cross section removed")
        elseif gname in GOLDEN_BOOLEAN_GUIDE && !isempty(aname)
            push!(deviations, "F4: $(dim)label $(repr(aname)) where Plots set the boolean $(repr(gname))")
        elseif isempty(expected) && !isempty(aname) && _has_contours(golden)
            push!(deviations, "F3: $(dim)label $(repr(aname)) kept where Plots dropped it")
        else
            push!(diffs, "$(at)$(dim)label $(repr(aname)) != $(repr(expected))")
        end
    end
    diffs
end

## ----------------------------------------------------------------------------------------------
## F4 and the grids bridge (step 8)
## ----------------------------------------------------------------------------------------------

# F4 — the electron/hole contribution plot loses its y guide under Plots. In
# `src/PlotRecipes/Waveform.jl`, `yguide --> u isa Unitful.Units{…} ? "Charge" : "Energy"` parses as
# `(yguide --> (u isa …)) ? … : …`, because an arrow operator binds tighter than a conditional. The
# recipe therefore writes the *boolean* into the guide and the golden spec records `"false"`. Naming
# the quantity is an accepted improvement (A1), the same way F3's kept labels are.
const GOLDEN_BOOLEAN_GUIDE = ("false", "true")

# F5 - a cross-section plot cannot have an axis in the dimension it cut away, so a golden guide that
# names the cut is Plots mislabelling the projection. `plot_electric_fieldlines` sets `yguide --> "z"`
# whatever the cut, which for a Cartesian `z =` cross section labels the y data "z"
# (`fieldlines_car_z`). Following the data is an accepted improvement (A1), like F3 and F4.
const TITLE_CUT = r"@\s*(?<dim>[^\s=]+)\s*="

_cut_dimension(title::AbstractString) =
    (m = match(TITLE_CUT, title); m === nothing ? "" : String(m[:dim]))

_is_grid(spec::CaseSpec) = spec.group == "grids"

"""
    _grid_panel_diffs(g, a, i; ext, deviations) -> Vector{String}

Bridge the three fields a grid panel moves (plan §2.7).

The Plots `Grid` recipe writes `"Grid point density - Axis N"` into the *x guide* of each subplot and
histograms the raw tick values, which carry no unit. So its golden spec has an empty title, that
string as the x label, and no axis unit at all, on an axis that is really a coordinate in internal
metres or radians. The Makie wrapper moves the string to the title and labels x with the coordinate
in display units (A1). What is checked instead:

* the title is exactly what Plots put on the x guide, panel for panel, so panels cannot reorder;
* the x label names a coordinate *and* a unit, and that unit is what
  [`_grid_extent_unit`](@ref) turns back into the internal unit the golden extent is measured in;
* the y label names the count Plots left blank.
"""
function _grid_panel_diffs(g, a, i::Int; ext, deviations)
    diffs = String[]
    g.xguide == a.title ||
        push!(diffs, "panel $i title $(repr(a.title)) != $(repr(g.xguide))")
    if _grid_extent_unit(a.xguide; ext) === nothing
        push!(diffs, "panel $i xlabel $(repr(a.xguide)) does not name a coordinate and a unit")
    else
        push!(deviations, "grids: panel $i labelled $(repr(a.xguide)) with " *
                          "$(repr(g.xguide)) moved to the title")
    end
    isempty(a.yguide) && push!(diffs, "panel $i ylabel is empty")
    diffs
end

# "r / mm" -> the *internal* unit those millimetres display, i.e. the unit the golden tick values are
# in. `display_unit` maps each internal unit to its display one, so this inverts it over the two
# internal units a grid axis can carry. `nothing` means the label named no unit at all.
const UNIT_LABEL = r"^(?<name>.+) / (?<unit>\S+)$"

function _grid_extent_unit(label::AbstractString; ext)
    m = match(UNIT_LABEL, label)
    m === nothing && return nothing
    shown = try
        _gate_unit(m[:unit])
    catch
        return nothing
    end
    for u in (ext.internal_length_unit, ext.internal_angle_unit)
        ext.display_unit(u) == shown && return u
    end
    nothing
end

"""
    _polar_data_diffs(golden_z, actual_z; ext, deviations) -> Vector{String}

Compare the data of a cylindrical `z` cross section, whose angular axis the two backends cannot
always hold identically.

Makie draws one flat quad per grid cell. A cell wider than `POLAR_MAX_CELL` is split into equal
parts, or the straight edges of its quad would cut across the disc as chords. That is a rendering
requirement, not a change to the data, so the angular row count may be a multiple of the golden
one. The radial extent and the values still have to match, and the block-mean digest too wherever
no cell was split.
"""
function _polar_data_diffs(gz, az; ext, deviations)
    diffs = String[]
    (gz === nothing || az === nothing) && return ["field z: absent on one side"]
    n_angles_golden, n_radial = gz.size[2], gz.size[1]
    az.size[2] == n_radial ||
        push!(diffs, "field z: radial size $(az.size[2]) != $(n_radial)")
    az.size[1] % n_angles_golden == 0 ||
        push!(diffs, "field z: angular size $(az.size[1]) is not a multiple of the grid's " *
                     "$(n_angles_golden)")
    PlotParitySpec._approx(gz.extrema, az.extrema; rtol = 1e-4) ||
        push!(diffs, "field z: extrema $(az.extrema) != $(gz.extrema)")
    if reverse(gz.digest_size) == az.digest_size
        # Loose, but it still catches a transposed or plain wrong slice, which is what the digest
        # is for.
        PlotParitySpec._approx(PlotParitySpec._transpose_digest(gz), az.digest; rtol = 0.05,
                atol = 0.05 * max(abs(gz.extrema[1]), abs(gz.extrema[2]), 1e-30)) ||
            push!(diffs, "field z: block-mean digest differs (orientation or values)")
    else
        push!(deviations, "polar: angular cells split into $(az.size[1]), " *
                          "digest not comparable")
    end
    diffs
end

# Equal data scaling, however the backend spells it: Plots' `aspect_ratio --> 1` records as `"1"` or
# `"1.0"` depending on whether the recipe wrote an `Int` or a `Float64`, an `Axis` carries
# `DataAspect()`, and an `Axis3` the symbol `:data`.
_equal_aspect(s::AbstractString) = s in ("1", "1.0", "data")

# Whether the golden plot puts values through a colour scale, and therefore had a colorbar next to
# it. `clims` is set by every scalar-field recipe, `colormap` is filled from the first
# gradient-coloured series, and the field series is the heatmap itself; a geometry plot has none.
_maps_colors(spec::CaseSpec) = spec.axes[1].clims !== nothing ||
    !isempty(spec.axes[1].colormap) || _field_series(spec) !== nothing

"""
    _extent_diffs(golden, actual; ext) -> Vector{String}

Compare the bounding box of everything drawn, dimension by dimension.

This is the only data check a geometry case gets. There is no heatmap to diff, and the two backends
hold a mesh completely differently: Plots emits one `:mesh3d` series per face carrying a `connections`
attribute, Makie one `mesh!` over all faces at once plus a `linesegments!`, so neither the series
count nor any single array is comparable. The *extent* is, and it is what catches a geometry drawn at
the wrong scale, in the wrong place, or with an axis swapped, all of which the title-and-label checks
above sail straight past.

It does not catch a wrong mesh of the right size. That is what the tier-2 gallery is for.
"""
function _extent_diffs(golden::CaseSpec, actual::CaseSpec; ext)
    diffs = String[]
    for panel in eachindex(golden.axes)
        g, a = golden.axes[panel], actual.axes[panel]
        at = length(golden.axes) == 1 ? "" : "panel $panel "
        # The unit the golden data is in. Plots recorded it on the axis, except for a grid, whose
        # ticks reached it as bare numbers; there it is recovered from the Makie label (§2.7).
        gunits = Any[g.xunit, g.yunit, g.zunit]
        _is_grid(golden) && (gunits[1] = _grid_extent_unit(a.xguide; ext))
        # A `PolarAxis` takes its angle in radians, the one place D7's display-unit convention does
        # not apply, so θ is compared exactly as the golden holds it.
        g.projection == "polar" && (gunits[1] = Unitful.NoUnits)
        for (i, dim) in enumerate(("x", "y", "z"))
            # Only a grid panel whose label named no unit, and `_grid_panel_diffs` has said so.
            gunits[i] === nothing && continue
            gb, ab = _data_bbox(golden, i, panel), _data_bbox(actual, i, panel)
            if gb === nothing || ab === nothing
                gb === ab || push!(diffs, "$(at)$dim extent: drawn on only one side " *
                                          "(golden: $(gb !== nothing), makie: $(ab !== nothing))")
                continue
            end
            expected = _expected_lims(gb, gunits[i]; ext)
            # A geometry is often flat in one dimension (an ellipse in z, a slice in its cut axis),
            # and a relative tolerance is meaningless against zero; the absolute one is scaled by the
            # size of the whole object so it stays a tight check on the dimensions that do have
            # extent.
            scale = maximum(abs, expected; init = 0.0)
            isapprox(ab, expected; rtol = 1e-4, atol = 1e-6 * max(scale, 1.0)) ||
                push!(diffs, "$(at)$dim extent $(ab) != $(expected)")
        end
    end
    diffs
end

# The interval the series of `spec` drawn into `panel` cover in dimension `i`, or `nothing` if
# nothing was drawn in it (a two-dimensional slice has no `z`).
function _data_bbox(spec::CaseSpec, i::Int, panel::Int = 1)
    lo, hi = Inf, -Inf
    for s in spec.series
        s.subplot == panel || continue
        d = (s.x, s.y, s.z)[i]
        (d === nothing || d.extrema === nothing) && continue
        lo, hi = min(lo, d.extrema[1]), max(hi, d.extrema[2])
    end
    isfinite(lo) ? Float64[lo, hi] : nothing
end

# A unit as the golden spec spells it. Some print under an abbreviation `uparse` cannot read back:
# the waveform signal unit `u"e_au"` prints as `"e"`, which is a symbol in no unit module.
# SolidStateDetectors' own config parser carries the alias table for exactly this, so the gate reuses
# it rather than inventing a second one.
function _gate_unit(s::AbstractString)
    isempty(s) && return Unitful.NoUnits
    get(() -> Unitful.uparse(s), SolidStateDetectors.unit_conversion, s)
end

# The cut position at the end of a title: `"… @ z = 0.0375 m"`, `"… @ φ = 0.0°"`.
const TITLE_POSITION = r"^(?<head>.*=\s)(?<value>-?[0-9.]+(?:e-?[0-9]+)?)\s*(?<unit>\S+)$"

"""
    _expected_title(golden_title; ext) -> String

The golden title with its cut position converted to display units.

Plots printed the position in the unit the caller passed — a bare `z = 0.04` meaning internal metres,
so `"@ z = 0.0375 m"` — while the Makie wrapper prints `display_unit` of it, `"@ z = 37.5 mm"`, so
that the title agrees with the axis labels. Everything else about the title (the field's name, the
cross-section symbol, the value itself) still has to match exactly; a title this cannot parse is
compared verbatim.
"""
function _expected_title(title::AbstractString; ext)
    m = match(TITLE_POSITION, title)
    m === nothing && return String(title)
    u = try
        _gate_unit(m[:unit])
    catch
        return String(title)
    end
    display = ext.display_unit(u)
    display == u && return String(title)
    q = parse(Float64, m[:value]) * u
    string(m[:head], round(display, uconvert(display, q), sigdigits = 3))
end

# The Plots guide plus the unit the Makie side displays it in: `("r", "m")` -> `"r / mm"`.
function _expected_label(name::AbstractString, unitstr::AbstractString; ext)
    isempty(name) && return ""
    u = ext.display_unit(_gate_unit(unitstr))
    u == Unitful.NoUnits ? String(name) : "$name / $u"
end

# Golden limits are in internal units; the Makie axis works in display units.
_expected_lims(lims, unitstr::AbstractString; ext) = _expected_lims(lims, _gate_unit(unitstr); ext)
_expected_lims(lims, u::Unitful.Units; ext) = Float64[ext.display_values(l, u) for l in lims]

# The series carrying the field itself. Plots always calls it a heatmap; on the Makie side a
# cylindrical `z` cross section is a mesh of polar cells instead, because a heatmap's axis-aligned
# cells do not survive the polar transform.
const FIELD_SERIES = ("heatmap", "surface", "polarcells")

_field_series(spec::CaseSpec) = findfirst(s -> s.seriestype in FIELD_SERIES, spec.series) |>
    i -> i === nothing ? nothing : spec.series[i]

# Plots names the overlay `contours`; on the Makie side `contour!` is itself a recipe, so the leaves
# under it are `lines`. Either way: a series that is not the field itself.
_has_contours(spec::CaseSpec) = any(s -> !(s.seriestype in FIELD_SERIES), spec.series)
