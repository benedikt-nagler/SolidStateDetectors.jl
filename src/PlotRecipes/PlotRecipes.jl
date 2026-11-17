# Backend-agnostic plotting helpers.
#
# The Plots.jl recipes that lived beside this file were removed in step 11 of the Makie migration.
# What remains is data preparation - cross sections, plot metadata, field-line sampling, legend
# labels - which the Makie extension in `ext/Makie/` consumes.
include("plot_utils.jl")
