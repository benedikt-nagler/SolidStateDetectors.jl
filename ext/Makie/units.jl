# Display units and axis labels.
#
# The two backends put units in opposite places. Plots' `unitformat = :slash` writes the unit into
# the *guide* (`"r / m"`) and leaves the tick labels bare; Makie's default `UnitfulConversion` does
# the reverse — it appends the unit to every tick label and leaves `ax.xlabel` empty.
#
# So each dimension gets a `UnitfulConversion(display_unit; units_in_label = false)` and the label is
# written by hand. Verified on Makie 0.24.13: data passed in internal units (metres, radians) with a
# `u"mm"` / `u"°"` conversion renders tick labels `0, 10, 20, 30, 40` under the label `"r / mm"` —
# exactly the GR look, and in a unit that reads better than GR's `0.035 m`.
#
# One trap worth spelling out: with a conversion installed, `ax.finallimits[]` and anything passed as
# `limits` are in **converted** (display) space, not in the data's own unit. `display_value` is the
# conversion to use when handing limits to an axis.

"""
    display_unit(u) -> Unitful.Units

The unit a quantity in internal unit `u` should be *shown* in.

Detector geometry is centimetre-scale, so lengths read far better in millimetres than in the internal
metres; angles read in degrees rather than radians; waveform times in nanoseconds rather than
seconds. Everything else (volts, V/m, V/m²) is already at a sensible scale and is left alone.
"""
function display_unit(u::Unitful.Units)
    u === Unitful.NoUnits && return u
    # Angles are dimensionless in Unitful, so they cannot be matched on dimension.
    (u == u"rad" || u == u"°") && return u"°"
    d = Unitful.dimension(u)
    d === Unitful.𝐋 && return u"mm"
    d === Unitful.𝐓 && return u"ns"
    u
end
display_unit(::Nothing) = Unitful.NoUnits

"""
    display_value(q, u = display_unit(unit(q))) -> Float64

A quantity as a bare number in its display unit — what an `Axis` with a `UnitfulConversion` expects
for `limits`, and what its `finallimits` are expressed in.
"""
display_value(q::Unitful.Quantity) = display_value(q, display_unit(Unitful.unit(q)))
display_value(q::Unitful.Quantity, u::Unitful.Units) = Float64(ustrip(uconvert(u, q)))
display_value(x::Real, ::Unitful.Units = Unitful.NoUnits) = Float64(x)

"""
    display_factor(u) -> Float64

The factor that takes a bare number in unit `u` to a bare number in `display_unit(u)` (1000 for
metres → millimetres, 180/π for radians → degrees, 1 for everything already at display scale).
"""
display_factor(u::Unitful.Units) = Float64(ustrip(uconvert(display_unit(u), 1.0 * u)))
display_factor(::Nothing) = 1.0

"""
    display_values(x, u) -> array

`x`, a bare number or array of bare numbers in unit `u`, rescaled to `display_unit(u)`.

Plot data has to be rescaled by hand rather than by an axis `UnitfulConversion`, because
**`heatmap` cannot take Unitful data on Makie 0.24**: `Heatmap`'s `convert_arguments` recurses
forever on a vector of `Quantity`s (`to_linspace` returns another `Quantity` range), and an axis
that carries a `UnitfulConversion` in turn rejects the bare numbers a heatmap does accept. So every
recipe here passes numbers already in display units and `ssd_axis` writes the unit into the label —
see the `dim_conversions` keyword there.
"""
function display_values(x::AbstractArray{<:Real}, u::Union{Unitful.Units, Nothing})
    f = display_factor(u)
    f == 1 && return collect(x)
    collect(x) .* convert(float(eltype(x)), f)   # `float`: a degree factor is not an integer
end
display_values(x::Real, u::Union{Unitful.Units, Nothing}) =
    x * convert(float(typeof(x)), display_factor(u))

"""
    unit_label(name, u; format = :slash) -> String

The axis or colorbar label for quantity `name` measured in `u`.

`:slash` gives `"r / mm"`, reproducing Plots' `unitformat = :slash`. `:parens` gives `"r (mm)"`,
which is the form Plots happens to use for *colorbar* titles. Unitless quantities get the bare name
in both cases.
"""
function unit_label(name::AbstractString, u::Unitful.Units; format::Symbol = :slash)
    (u == Unitful.NoUnits || isempty(string(u))) && return String(name)
    format === :slash && return "$name / $u"
    format === :parens && return "$name ($u)"
    throw(ArgumentError("unknown unit label format $(repr(format)); use :slash or :parens"))
end
unit_label(name::AbstractString, ::Nothing; kwargs...) = String(name)

"""
    unit_conversion(u) -> Makie.UnitfulConversion or nothing

The `dim*_conversion` for an axis whose data carries unit `u`: it rescales the data to
`display_unit(u)` and keeps the unit *out* of the tick labels, because `unit_label` puts it in the
axis label instead. Unitless dimensions get `nothing`, which is Makie's own default for the
attribute (`Makie.automatic` is *not* — an `Axis` refuses to convert it to an `AbstractDimConversion`).
"""
function unit_conversion(u::Unitful.Units)
    u == Unitful.NoUnits && return nothing
    Makie.UnitfulConversion(display_unit(u); units_in_label = false)
end
unit_conversion(::Nothing) = nothing
