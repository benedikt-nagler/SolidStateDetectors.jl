# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

module SolidStateDetectorsCairoMakieExt

# CairoMakie paints a mesh one triangle at a time, each as its own Cairo mesh pattern. Neighbouring
# patterns overlap by a fraction of a pixel along their shared edge, which is invisible on an opaque
# mesh and draws every such edge as a darker line on a translucent one. Detector geometry is
# translucent by design (a contact at 20 % opacity so the crystal shows through), so every `:csg` and
# `:mesh3d` picture came out covered in the triangulation of its faces.
#
# The faces of a geometry plot are flat-coloured and unshaded, so they do not need a mesh pattern at
# all. This extension fills each triangle as a plain path with antialiasing off, which makes
# neighbours share their edge pixel-exactly: the surface comes out as one even tone.
#
# It hooks in through `CairoMakie.draw_plot`, the same function CairoMakie's own `poly` and `band`
# overrides use, and reads a few CairoMakie internals to project the vertices. Should any of that
# change, `draw_flat_mesh` returns `false` and the mesh is drawn the stock way: the seams come back,
# nothing breaks.

using Makie
import CairoMakie
import CairoMakie: Cairo
import SolidStateDetectors

const CSGPlot = Makie.Plot{SolidStateDetectors.csgplot}
const FieldPlot = Union{Makie.Plot{SolidStateDetectors.scalarpotentialplot},
                        Makie.Plot{SolidStateDetectors.electricfieldplot},
                        Makie.Plot{SolidStateDetectors.fieldlinesplot}}

# CairoMakie draws the plots it calls atomic and recurses into every other one. Declaring a geometry
# plot atomic is what routes it through `draw_plot` below, as CairoMakie does for its own `poly`.
CairoMakie.is_cairomakie_atomic_plot(::CSGPlot) = true

# A field plot carries a `rasterize` attribute, and CairoMakie takes any plot that has one for an
# atomic it should know how to draw. So it does a plot with no children, which is what field lines
# on a cross section no contact meets come to. It warns on every render and then draws the children
# anyway. There is nothing to draw for the composite itself.
CairoMakie.draw_atomic(::Makie.Scene, ::CairoMakie.Screen, ::FieldPlot) = nothing

function CairoMakie.draw_plot(scene::Makie.Scene, screen::CairoMakie.Screen, p::CSGPlot)
    Makie.to_value(get(p, :visible, true)) || return
    children = p.plots
    for idx in sortperm(Makie.zvalue2d.(children))
        child = children[idx]
        (child isa Makie.Mesh && draw_flat_mesh(scene, screen, child)) ||
            CairoMakie.draw_plot(scene, screen, child)
    end
    return
end

"""
    draw_flat_mesh(scene, screen, plot::Makie.Mesh) -> Bool

Fill the triangles of an unshaded, single-colour 3D mesh back to front, without seams.

Returns `false`, having drawn nothing, for any mesh this does not cover or when the CairoMakie
internals it relies on are not there; the caller then draws the mesh the stock way.
"""
function draw_flat_mesh(scene, screen, plot)
    local pts, faces, color, order
    try
        Makie.to_value(get(plot, :visible, true)) || return true
        plot.shading[] && return false
        color = CairoMakie.compute_colors(plot)
        color isa Makie.Colorant || return false
        pts = CairoMakie.cairo_project_to_screen(plot; output_type = Makie.Point3f)
        (pts isa AbstractVector && eltype(pts) <: Makie.Point3) || return false
        faces = plot.faces[]
        keep = _unclipped_faces(scene, plot, faces)
        depth(f) = (pts[f[1]][3] + pts[f[2]][3] + pts[f[3]][3]) / 3
        order = sort!(keep; by = i -> depth(faces[i]), rev = true)
    catch
        return false
    end

    ctx = screen.context
    Cairo.save(ctx)
    try
        Cairo.set_antialias(ctx, Cairo.ANTIALIAS_NONE)
        Cairo.set_source_rgba(ctx, Makie.red(color), Makie.green(color), Makie.blue(color),
                              Makie.alpha(color))
        for i in order
            f = faces[i]
            a, b, c = pts[f[1]], pts[f[2]], pts[f[3]]
            (all(isfinite, a) && all(isfinite, b) && all(isfinite, c)) || continue
            Cairo.move_to(ctx, a[1], a[2])
            Cairo.line_to(ctx, b[1], b[2])
            Cairo.line_to(ctx, c[1], c[2])
            Cairo.close_path(ctx)
            Cairo.fill(ctx)
        end
    finally
        Cairo.restore(ctx)
    end
    return true
end

# The faces left after clipping, as CairoMakie decides it: a face with a vertex outside the clip
# planes of its axis is dropped whole.
function _unclipped_faces(scene, plot, faces)
    planes = plot.clip_planes[]
    (isempty(planes) || !Makie.is_data_space(plot.space[])) && return collect(eachindex(faces))
    world = Makie.apply_model(plot.model_f32c[], plot.positions_transformed_f32c[])
    valid = Bool[CairoMakie.is_visible(planes, pt) for pt in world]
    Int[i for i in eachindex(faces) if all(v -> valid[v], faces[i])]
end

end # module
