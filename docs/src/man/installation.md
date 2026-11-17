# Installation

This package is a registered package.

Install via

```julia
using Pkg; pkg"add SolidStateDetectors"
```

## Visualization / Plotting (Optional)

This package provides plotting functions for its geometries and simulation results, built on
[Makie.jl](https://docs.makie.org). They live in a package extension, so they become available as
soon as a Makie backend is loaded.

Install one via

```julia
using Pkg; pkg"add CairoMakie"
```

Load a [Makie.jl](https://docs.makie.org) backend via

```julia
using CairoMakie
```

CairoMakie writes vector graphics and needs no display, which suits documents and headless machines.
[GLMakie.jl](https://docs.makie.org/stable/explanations/backends/glmakie) opens an interactive window
instead and renders three-dimensional detector geometry better; either one activates the plotting
functions.

For more information about the plotting functions of this package look up the [Plotting](@ref) section.

This documentation was built with
```@example
using CairoMakie # hide
print("CairoMakie: v$(pkgversion(CairoMakie)) - Makie: v$(pkgversion(CairoMakie.Makie))") # hide
```

## GPU Support in Field Calculations

The [Electric Potential](@ref) and individual [Weighting Potentials](@ref) can also be calculated on GPUs. SolidStateDetectors.jl uses [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) for GPU support.

In order to use your GPU, the Julia Packages [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) or [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl) have to be installed and loaded. SolidStateDetectors.jl may run on other GPU types supported by KernelAbstractions.jl, but has not been tested with [oneAPI.jl](https://github.com/JuliaGPU/oneAPI.jl) or [Metal.jl](https://github.com/JuliaGPU/Metal.jl) yet.

For NVIDIA GPUs use:

```julia
using CUDA, SolidStateDetectors
# *Only* for KernelAbstractions < v0.9:
# using CUDAKernels
```

For AMD GPUs use:

```julia
using AMDGPU, SolidStateDetectors
# *Only* for KernelAbstractions < v0.9:
# using ROCKernels
```

Then, in any field calculation ([`calculate_electric_potential!`](@ref), [`calculate_weighting_potential!`](@ref), [`simulate!(::Simulation)`](@ref)) the keyword `device_array_type` can be set to choose the device on which the calculations should be performed.
The possibilities are:
```julia
device_array_type = Array # -> CPU (default)
device_array_type = CuArray # -> NVIDIA GPU
device_array_type = ROCArray # -> AMD GPU
```

### Example (NVIDIA)

```julia
using CUDA, SolidStateDetectors
# *Only* for KernelAbstractions < v0.9:
# using CUDAKernels

sim = Simulation(SSD_examples[:CGD])
calculate_electric_potential!( 
    sim, 
    device_array_type = CuArray, 
    refinement_limits = [0.2, 0.1, 0.05],
    depletion_handling = true
)
```

!!! note
    The AMD backend was not yet tested due to lack of an AMD GPU (we are working on that).
