
"""
Structural AbstractPathmodel. 
Consists of a graph from StenoGraphs structural (structural model). 

## Arguments

- `structural_model`: Graph from StenoGraphs package. Defines the latent relations between the factors of measurement_model.  

```julia
using StenoGraphs
using Taxonomy

graph = @StenoGraph begin
    # latent regressions
    fac1 → fac2
end

Structural(structural_model = graph)

# output

Structural
   n_sample: Judgement{Missing}
   structural_model: Judgement{Vector{DirectedEdge{SimpleNode{Symbol}, SimpleNode{Symbol}}}}

```
"""
struct Structural <: AbstractPathmodel
    n_sample::Judgement{<: Int}
    structural_model::Judgement{ <: Union{<:AbstractArray{<: StenoGraphs.AbstractEdge}, Missing}}
    Structural(n_sample, structural_model) = new(J(n_sample), J(structural_model))
end

function Structural(;n_sample = missing,
    structural_model) # hier basic CFA model festlegen
    Structural(n_sample, structural_model)
end


function ManifestPathmodel(;n_sample, kwargs...)
    Structural(;n_sample = n_sample, kwargs...)
end
