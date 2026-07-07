# Backend-agnostic plotting helpers.
#
# These functions contain no RecipesBase/Plots (or Makie) code: they only do data preparation
# (cross-section lookup, plot metadata, field-line sampling) that is shared by both the existing
# Plots.jl recipes and the Makie extension. They live here (rather than inside the `@recipe` files)
# so they survive the eventual removal of the Plots.jl recipes.

_get_potential_plot_information(sp::ScalarPotential, ::Bool) = _get_potential_plot_information(sp) # ignore point type bool for all other potentials

function _get_potential_plot_information(p::ElectricPotential)::Tuple{Symbol, Tuple{Quantity, Quantity}, String, Unitful.Units}
    return (:viridis, extrema(p.data).*internal_voltage_unit, "Electric Potential", internal_voltage_unit)
end

function _get_potential_plot_information(::WeightingPotential{T})::Tuple{Symbol, Tuple{T,T}, String, Unitful.Units} where {T}
    return (:viridis, (zero(T), one(T)), "Weighting Potential", Unitful.NoUnits)
end

function _get_potential_plot_information(ecd::EffectiveChargeDensity)::Tuple{Symbol, Tuple{Quantity, Quantity}, String, Unitful.Units}
    return (:inferno, extrema(ecd.data).*u"V/m^2", "Effective Charge Density", u"V/m^2")
end

function _get_potential_plot_information(::PointTypes{T}, all_point_types::Bool)::Tuple{Symbol, Tuple{T,T}, String, Unitful.Units} where {T}
    return (:viridis, (zero(0), (all_point_types ? 31 : 7)*one(T)), "Point Type Map", Unitful.NoUnits)
end

function _get_potential_plot_information(::ImpurityScale{T})::Tuple{Symbol, Tuple{T,T}, String, Unitful.Units} where {T}
    return (:viridis, (zero(T), one(T)), "Impurity Scale", Unitful.NoUnits)
end

function get_crosssection_idx_and_value(grid::Grid, ::Any, ::Any, ::Any)
    error(ArgumentError, ": The cross section your are trying to plot is not given in the right format.\nPlease pass a number with matching or no units.")
end


## Cross sections of cylindrical grids

function get_crosssection_idx_and_value(grid::CylindricalGrid{T}, r::Union{Real, LengthQuantity, Missing}, φ::Union{Real, AngleQuantity, Missing}, z::Union{Real, LengthQuantity, Missing})::Tuple{Symbol,Int,T,Unitful.Units} where {T <: SSDFloat}

    cross_section::Symbol, idx::Int, units::Unitful.Units = if ismissing(φ) && ismissing(r) && ismissing(z)
        return get_crosssection_idx_and_value(grid, r, T(0.0), z)
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        φ_rad::T = T(to_internal_units(φ isa Real ? deg2rad(φ) : φ))
        while !(grid.φ.interval.left <= φ_rad <= grid.φ.interval.right) && grid.φ.interval.right != grid.φ.interval.left
            if φ_rad > grid.φ.interval.right
                φ_rad -= width(grid.φ.interval)
            elseif φ_rad < grid.φ.interval.left
                φ_rad += width(grid.φ.interval)
            end
        end
        :φ, searchsortednearest(grid.φ, φ_rad), φ isa Real ? u"°" : unit(φ)
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(grid.r, T(to_internal_units(r))), r isa Real ? internal_length_unit : unit(r)
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(grid.z, T(to_internal_units(z))), z isa Real ? internal_length_unit : unit(z)
    else
        error(ArgumentError, ": Only one of the keywords `r, φ, z` is allowed.")
    end

    value::T = if cross_section == :φ
        rad2deg(grid.φ[idx])
    elseif cross_section == :r
        grid.r[idx]
    elseif cross_section == :z
        grid.z[idx]
    end

    cross_section, idx, value, units
end


## Cross sections of Cartesian grids

function get_crosssection_idx_and_value(grid::CartesianGrid3D{T}, x::Union{Real, LengthQuantity, Missing}, y::Union{Real, LengthQuantity, Missing}, z::Union{Real, LengthQuantity, Missing})::Tuple{Symbol,Int,T,Unitful.Units} where {T <: SSDFloat}

    cross_section::Symbol, idx::Int, units::Unitful.Units = if ismissing(x) && ismissing(y) && ismissing(z)
        return get_crosssection_idx_and_value(grid, T(0.0), y, z)
    elseif !ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(grid.x, T(to_internal_units(x))), x isa Real ? internal_length_unit : unit(x)
    elseif ismissing(x) && !ismissing(y) && ismissing(z)
        :y, searchsortednearest(grid.y, T(to_internal_units(y))), y isa Real ? internal_length_unit : unit(y)
    elseif ismissing(x) && ismissing(y) && !ismissing(z)
        :z, searchsortednearest(grid.z, T(to_internal_units(z))), z isa Real ? internal_length_unit : unit(z)
    else
        error(ArgumentError, ": Only one of the keywords `x, y, z` is allowed.")
    end
    value::T = grid[cross_section][idx]
    cross_section, idx, value, units
end


## Field-line sampling lines (used by electric field line plotting)

function get_sample_lines(dim_symbol::Symbol, v::T, grid::CartesianGrid3D{T}, sampling::T)::Vector{ConstructiveSolidGeometry.Line} where {T}

    xrange = range(round.(endpoints(grid.x.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    yrange = range(round.(endpoints(grid.y.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    zrange = range(round.(endpoints(grid.z.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)

    return if dim_symbol == :x
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(v,y,0), CartesianVector{T}(0,0,1)) for y in yrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(v,0,z), CartesianVector{T}(0,1,0)) for z in zrange]
        )
    elseif dim_symbol == :y
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(x,v,0), CartesianVector{T}(0,0,1)) for x in xrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,v,z), CartesianVector{T}(1,0,0)) for z in zrange]
        )

    else # dim_symbol == :z
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(x,0,v), CartesianVector{T}(0,1,0)) for x in xrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,y,v), CartesianVector{T}(1,0,0)) for y in yrange]
        )
    end
end


function get_sample_lines(dim_symbol::Symbol, v::T, grid::CylindricalGrid{T}, sampling::T)::Vector{ConstructiveSolidGeometry.Line} where {T}

    φsampling = 2π * sampling / (dim_symbol == :r ? v : grid.r.interval.right) # or use sampling together with the unit ?
    rrange = range(round.((-grid.r.interval.right,grid.r.interval.right)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    φrange = range(round.(endpoints(grid.φ.interval)./φsampling, (RoundDown, RoundUp)).*φsampling..., step = φsampling)
    zrange = range(round.(endpoints(grid.z.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)

    return if dim_symbol == :φ
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(r*cos(v),r*sin(v),0), CartesianVector{T}(0,0,1)) for r in rrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,0,z), CartesianVector{T}(cos(v),sin(v),0)) for z in zrange]
        )
    elseif dim_symbol == :z
        [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,0,v), CartesianVector{T}(cos(φ), sin(φ),0)) for φ in φrange]
    else # dim_symbol == :r
        [ConstructiveSolidGeometry.Line(CartesianPoint{T}(v*cos(φ),v*sin(φ),0), CartesianVector{T}(0,0,1)) for φ in φrange]
    end
end
