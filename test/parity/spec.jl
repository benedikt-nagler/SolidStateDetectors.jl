# Plot-parity specification: the backend-neutral description of "what a plot is".
#
# Part of the Plots.jl -> Makie.jl migration harness. Two producers
# write these specs:
#
#   * `devscripts/extract_plots_spec.jl` introspects a Plots.jl `Plot` object (tier 0) and writes one
#     JSON per case into `test/parity/golden/`. Those JSONs are committed and stay valid after the
#     Plots.jl recipes are deleted — they *are* the frozen specification of current behaviour.
#   * the tier-1 assertions in `test/test_plotting.jl` introspect the Makie figure of the same case
#     and compare it against the committed golden spec.
#
# This file therefore must not depend on either backend: it only defines the data model, the JSON
# round trip and the comparison helpers.
#
# Numeric payloads are not stored verbatim (that would be megabytes of JSON per run). Each array is
# reduced to a `DataSpec`: shape, unit, extrema, a content hash, and a small block-mean `digest`.
# The digest is what makes a *transposed* array detectable — the single most likely silent error of
# this migration (plan §1.7) — because block means transpose with the array while shape, extrema and
# even the hash of the flattened data may not change.

module PlotParitySpec

using JSON
using Statistics: mean
using Unitful

export DataSpec, SeriesSpec, AxisSpec, CaseSpec
export data_spec, write_spec, read_spec, spec_path, has_spec, compare_data, compare_colormap

const SPEC_VERSION = 1

# Number of blocks per dimension the block-mean digest is reduced to.
const DIGEST_BLOCKS = 8

## ----------------------------------------------------------------------------------------------
## Data model
## ----------------------------------------------------------------------------------------------

"""
    DataSpec

Fingerprint of one numeric plot argument (a series' `x`, `y` or `z`).

* `size`: shape of the array (`[n]` for a vector).
* `unit`: the Unitful unit as a string, `""` if the data is unitless.
* `extrema`: `[min, max]` of the unit-stripped values, `nothing` if the array is empty or all-missing.
* `digest`: block means over a grid of at most `DIGEST_BLOCKS` blocks per dimension, flattened in
  column-major order. Transposing the array transposes the digest, so comparing digests catches an
  orientation error that shape and extrema would not.
* `digest_size`: shape of the block grid the digest was flattened from.
* `hash`: hash of the rounded values; equal hashes mean identical data, unequal hashes mean nothing
  in particular (the two backends reach the same plot through different arithmetic).
* `n_missing`: number of `missing`/`NaN` entries, which waveform and drift-path data legitimately have.
"""
struct DataSpec
    size::Vector{Int}
    unit::String
    extrema::Union{Nothing, Vector{Float64}}
    digest::Vector{Float64}
    digest_size::Vector{Int}
    hash::String
    n_missing::Int
end

"""
    SeriesSpec

One drawn series. In Plots this is an element of `plt.series_list`; in Makie it is one child plot of
the recipe's plot tree.
"""
struct SeriesSpec
    subplot::Int
    seriestype::String
    label::String
    color::String
    linecolor::String
    linewidth::Union{Nothing, Float64}
    alpha::Union{Nothing, Float64}
    fillalpha::Union{Nothing, Float64}
    markershape::String
    markersize::Union{Nothing, Float64}
    n_connections::Int
    x::Union{Nothing, DataSpec}
    y::Union{Nothing, DataSpec}
    z::Union{Nothing, DataSpec}
end

"""
    AxisSpec

Everything a Plots recipe sets at subplot level and a Makie recipe *cannot* set at all (plan §1.1):
title, guides, limits, colour scale, colorbar and legend. This is the part of the spec that only the
figure-level wrapper functions can satisfy, and the reason smoke tests are banned as a gate.

Note where the units are: Plots' Unitful integration strips units off the series data and records
them per axis (`sp[:xaxis][:unit]`), leaving the guide bare (`"r"`, with `unitformat = :slash`
appending `/ m` only at render time). So `xunit`/`yunit`/`zunit` hold the unit and the `DataSpec`s of
the series are plain numbers in that unit. `xlims`/`ylims`/`clims` are likewise unit-stripped.
"""
struct AxisSpec
    title::String
    xguide::String
    yguide::String
    zguide::String
    xunit::String
    yunit::String
    zunit::String
    xlims::Union{Nothing, Vector{Float64}}
    ylims::Union{Nothing, Vector{Float64}}
    clims::Union{Nothing, Vector{Float64}}
    colorbar::String
    colorbar_title::String
    colormap::String
    aspect_ratio::String
    projection::String
    legend_position::String
    legend_entries::Vector{String}
end

"""
    CaseSpec

The full spec of one case of the matrix (`test/parity/cases.jl`).

* `error` is non-empty when the case could not be plotted at all — recorded rather than thrown, so
  one broken case does not stop a run.
* `note` is non-empty when the recorded behaviour is known *not* to be a valid parity reference (see
  `PlotParityCases.NOTE_R_ALIAS`). Tier-1 must not diff against a spec that carries a note.
"""
struct CaseSpec
    version::Int
    id::String
    group::String
    call::String
    backend::String
    axes::Vector{AxisSpec}
    series::Vector{SeriesSpec}
    error::String
    note::String
end

CaseSpec(; id, group, call, backend, axes = AxisSpec[], series = SeriesSpec[], error = "", note = "") =
    CaseSpec(SPEC_VERSION, id, group, call, backend, axes, series, error, note)

## ----------------------------------------------------------------------------------------------
## Building a DataSpec from an array
## ----------------------------------------------------------------------------------------------

_unit_string(::Any) = ""
_unit_string(u::Unitful.Units) = string(u)

# Unwrap whatever a backend hands us into (plain Float64 array, unit string).
function _strip(a::AbstractArray)
    el = eltype(a)
    if el <: Unitful.Quantity
        u = Unitful.unit(el)
        return Float64.(ustrip.(u, a)), _unit_string(u)
    elseif el <: Union{Missing, <:Unitful.Quantity}
        u = Unitful.unit(nonmissingtype(el))
        return [ismissing(x) ? NaN : Float64(ustrip(u, x)) for x in a], _unit_string(u)
    elseif el <: Real || el <: Union{Missing, <:Real}
        return [ismissing(x) ? NaN : Float64(x) for x in a], ""
    elseif el <: Integer
        return Float64.(a), ""
    else
        # PointType (a UInt8 wrapper), Bool, and anything else convertible one element at a time.
        return [Float64(x) for x in a], ""
    end
end

"""
    data_spec(a) -> DataSpec
    data_spec(nothing) -> nothing

Fingerprint a plot argument. Returns `nothing` for `nothing` and for anything not reducible to
numbers (a colour vector, say), so a caller can record "this series has no such argument".
"""
data_spec(::Nothing) = nothing
data_spec(a) = data_spec(collect(a))
function data_spec(a::AbstractArray)
    vals, unitstr = try
        _strip(a)
    catch
        return nothing
    end
    vals isa AbstractArray{Float64} || return nothing
    n_missing = count(isnan, vals)
    finite = filter(isfinite, vec(vals))
    ex = isempty(finite) ? nothing : Float64[minimum(finite), maximum(finite)]
    dig, digsize = _digest(vals)
    h = string(hash(round.(vec(vals), sigdigits = 10)), base = 16)
    DataSpec(collect(Int, size(vals)), unitstr, ex, dig, collect(Int, digsize), h, n_missing)
end
data_spec(x::Number) = data_spec([x])

# Block-mean reduction to at most DIGEST_BLOCKS blocks per dimension. NaNs are skipped; an all-NaN
# block becomes NaN, which round-trips through JSON as `null` and compares equal to another NaN
# block (see `compare_data`).
function _digest(a::AbstractArray{Float64, N}) where {N}
    dims = size(a)
    nblocks = ntuple(i -> min(DIGEST_BLOCKS, dims[i]), N)
    any(iszero, dims) && return Float64[], collect(Int, nblocks)
    edges = ntuple(N) do i
        [round(Int, (b - 1) * dims[i] / nblocks[i]) + 1 for b in 1:(nblocks[i] + 1)]
    end
    edges = ntuple(i -> vcat(edges[i][1:nblocks[i]], dims[i] + 1), N)
    out = Vector{Float64}(undef, prod(nblocks))
    for (k, idx) in enumerate(CartesianIndices(nblocks))
        ranges = ntuple(i -> edges[i][idx[i]]:(edges[i][idx[i] + 1] - 1), N)
        block = @view a[ranges...]
        vals = filter(isfinite, vec(block))
        out[k] = isempty(vals) ? NaN : mean(vals)
    end
    out, collect(Int, nblocks)
end
# Anything else is converted once and rejected if that did not produce plain numbers — a vector of
# `Point{3, Float32}`, say, survives `Float64.(…)` as a vector of `Point{3, Float64}`, and a
# recursive fallback would then loop until the stack overflows.
function _digest(a::AbstractArray)
    b = Float64.(a)
    b isa AbstractArray{Float64} ||
        throw(ArgumentError("cannot digest an array with eltype $(eltype(a))"))
    _digest(b)
end

## ----------------------------------------------------------------------------------------------
## Comparison
## ----------------------------------------------------------------------------------------------

_approx(a::Nothing, b::Nothing; kw...) = true
_approx(a, b::Nothing; kw...) = false
_approx(a::Nothing, b; kw...) = false
_approx(a::Real, b::Real; rtol = 1e-6, atol = 0.0) =
    (isnan(a) && isnan(b)) || isapprox(a, b; rtol, atol)
_approx(a::AbstractVector, b::AbstractVector; kw...) =
    length(a) == length(b) && all(_approx(x, y; kw...) for (x, y) in zip(a, b))

"""
    compare_data(golden, actual; rtol = 1e-4, transposed = false) -> Vector{String}

Diff two `DataSpec`s, returning one human-readable line per mismatch (empty means they match).

`transposed = true` compares `actual` against `golden` **transposed** — which is how a Makie heatmap
`z` must be compared against the Plots golden spec, since Plots wants `(length(y), length(x))` and
Makie wants `(length(x), length(y))` (plan §1.7). Comparison is on shape, unit, extrema and the
block-mean digest; the hash is informational only, because the two backends reach the same picture
through different arithmetic.
"""
function compare_data(golden::Union{Nothing, DataSpec}, actual::Union{Nothing, DataSpec};
                      rtol::Real = 1e-4, transposed::Bool = false, label::AbstractString = "data")
    diffs = String[]
    if golden === nothing || actual === nothing
        golden === actual || push!(diffs, "$label: one side is absent (golden=$(golden === nothing ? "absent" : "present"), actual=$(actual === nothing ? "absent" : "present"))")
        return diffs
    end
    gsize, gdigest, gdigsize = if transposed
        reverse(golden.size), _transpose_digest(golden), reverse(golden.digest_size)
    else
        golden.size, golden.digest, golden.digest_size
    end
    gsize == actual.size || push!(diffs, "$label: size $(actual.size) != $(gsize)$(transposed ? " (golden transposed)" : "")")
    golden.unit == actual.unit || push!(diffs, "$label: unit \"$(actual.unit)\" != \"$(golden.unit)\"")
    _approx(golden.extrema, actual.extrema; rtol) ||
        push!(diffs, "$label: extrema $(actual.extrema) != $(golden.extrema)")
    if gdigsize == actual.digest_size
        _approx(gdigest, actual.digest; rtol, atol = rtol * _scale(golden)) ||
            push!(diffs, "$label: block-mean digest differs (orientation or values)")
    else
        push!(diffs, "$label: digest shape $(actual.digest_size) != $(gdigsize)")
    end
    diffs
end

_scale(d::DataSpec) = d.extrema === nothing ? 1.0 : max(abs(d.extrema[1]), abs(d.extrema[2]), 1e-30)

"""
    compare_colormap(golden, actual; atol = 12/255) -> Vector{String}

Diff two colormap strings. Both backends record a colormap as three sampled stops
(`"grad(#440154,#21908c,#fde725)"`), but they sample the underlying ColorScheme differently — Plots
interpolates the gradient, Makie indexes a 256-entry vector — so the midpoint stop can differ by a
step or two. Hence a tolerance rather than string equality; identical strings short-circuit.
"""
function compare_colormap(golden::AbstractString, actual::AbstractString; atol::Real = 12 / 255)
    golden == actual && return String[]
    g, a = _parse_stops(golden), _parse_stops(actual)
    (g === nothing || a === nothing || length(g) != length(a)) &&
        return ["colormap: \"$actual\" != \"$golden\""]
    maxdiff = maximum(maximum(abs.(x .- y)) for (x, y) in zip(g, a))
    maxdiff <= atol ? String[] :
        ["colormap: \"$actual\" != \"$golden\" (max channel difference $(round(maxdiff, digits = 3)))"]
end

function _parse_stops(s::AbstractString)
    m = match(r"^grad\((.*)\)$", s)
    m === nothing && return nothing
    stops = Vector{Float64}[]
    for hex in split(m.captures[1], ",")
        h = strip(hex)
        (startswith(h, "#") && length(h) in (7, 9)) || return nothing
        push!(stops, [parse(Int, h[i:(i + 1)], base = 16) / 255 for i in 2:2:6])
    end
    stops
end

# Reinterpret a flattened 2D digest as its transpose. 1D digests are unaffected.
function _transpose_digest(d::DataSpec)
    length(d.digest_size) == 2 || return d.digest
    m = reshape(d.digest, d.digest_size[1], d.digest_size[2])
    vec(permutedims(m))
end

## ----------------------------------------------------------------------------------------------
## JSON round trip
## ----------------------------------------------------------------------------------------------

const GOLDEN_DIR = joinpath(@__DIR__, "golden")

"""
    spec_path(id; dir = GOLDEN_DIR)

Path of the JSON file holding the spec of case `id`.
"""
spec_path(id::AbstractString; dir::AbstractString = GOLDEN_DIR) = joinpath(dir, "$(id).json")

_json(x::Nothing) = nothing
# NaN/Inf are not valid JSON numbers. Values are rounded to 8 significant digits: far tighter than
# the 1e-4 comparison tolerance, but it keeps the committed files a third smaller and stops the last
# float digit from churning the diff between runs.
_json(x::Real) = isfinite(x) ? round(Float64(x), sigdigits = 8) : nothing
_json(x::Integer) = x
_json(x::AbstractString) = String(x)
_json(x::AbstractVector) = Any[_json(v) for v in x]
_json(x::T) where {T} = Dict{String, Any}(string(f) => _json(getfield(x, f)) for f in fieldnames(T))

_num(v, default = NaN) = v === nothing ? default : Float64(v)
_numvec(v) = v === nothing ? nothing : Float64[_num(x) for x in v]
_str(v) = v === nothing ? "" : String(v)

_dataspec(::Nothing) = nothing
_dataspec(d::AbstractDict) = DataSpec(
    Int[Int(x) for x in d["size"]], _str(d["unit"]), _numvec(d["extrema"]),
    Float64[_num(x) for x in d["digest"]], Int[Int(x) for x in d["digest_size"]],
    _str(d["hash"]), Int(d["n_missing"]))

_optnum(v) = v === nothing ? nothing : Float64(v)

_seriesspec(d::AbstractDict) = SeriesSpec(
    Int(d["subplot"]), _str(d["seriestype"]), _str(d["label"]), _str(d["color"]), _str(d["linecolor"]),
    _optnum(d["linewidth"]), _optnum(d["alpha"]), _optnum(d["fillalpha"]),
    _str(d["markershape"]), _optnum(d["markersize"]), Int(d["n_connections"]),
    _dataspec(d["x"]), _dataspec(d["y"]), _dataspec(d["z"]))

_axisspec(d::AbstractDict) = AxisSpec(
    _str(d["title"]), _str(d["xguide"]), _str(d["yguide"]), _str(d["zguide"]),
    _str(d["xunit"]), _str(d["yunit"]), _str(d["zunit"]),
    _numvec(d["xlims"]), _numvec(d["ylims"]), _numvec(d["clims"]),
    _str(d["colorbar"]), _str(d["colorbar_title"]), _str(d["colormap"]),
    _str(d["aspect_ratio"]), _str(d["projection"]),
    _str(d["legend_position"]), String[_str(x) for x in d["legend_entries"]])

"""
    write_spec(spec; dir = GOLDEN_DIR) -> String

Write `spec` as pretty-printed JSON (stable key order, so the files diff cleanly) and return the path.
"""
function write_spec(spec::CaseSpec; dir::AbstractString = GOLDEN_DIR)
    mkpath(dir)
    path = spec_path(spec.id; dir)
    open(path, "w") do io
        JSON.print(io, _json(spec), 2)
    end
    path
end

"""
    read_spec(id; dir = GOLDEN_DIR) -> CaseSpec

Read the committed golden spec of case `id`.
"""
function read_spec(id::AbstractString; dir::AbstractString = GOLDEN_DIR)
    d = JSON.parsefile(spec_path(id; dir))
    CaseSpec(Int(d["version"]), _str(d["id"]), _str(d["group"]), _str(d["call"]), _str(d["backend"]),
             AxisSpec[_axisspec(a) for a in d["axes"]],
             SeriesSpec[_seriesspec(s) for s in d["series"]],
             _str(d["error"]), _str(get(d, "note", "")))
end

"""
    has_spec(id; dir = GOLDEN_DIR)

Whether a golden spec exists for case `id`.
"""
has_spec(id::AbstractString; dir::AbstractString = GOLDEN_DIR) = isfile(spec_path(id; dir))

end # module
