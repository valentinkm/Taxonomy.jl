"""
CFA Factor Taxon. 
Building Block for CFA Taxonomy. Multiple Factors can be combined to a CFA. 

## Arguments

- `n_sample`: Number of observed cases. May be between different taxons from same paper sometime, e.g. multigroup models.
- `n_variables`: Number of variables (possibly observed/manifest).
- `loadings`: Vector of loadings, one for each item. 
- `factor_variance`: Variance of the factor.
- `error_variances`: Vector of variances of the respective errors
- `error_covariances_within`: Vector of covariances within factor.
- `error_covariances_between`: Vector of covariances the factor shares with a different factor. 
- `crossloadings_incoming`: Vector of crossloadings coming from other factors. They should be lower than the loading coming to the item from this factor.  
- `crossloadings_outgoing`: Vector of crossloadings going to other items which have higher loadings from other factors. 

```jldoctest
CFAFactor(n_variables = 2, loadings = [1, 0.4], factor_variance = 0.6)

# output

CFAFactor
   n_sample: Judgement{Missing}
   n_variables: Judgement{Int64}
   loadings: Judgement{Vector{Float64}}
   factor_variance: Judgement{Float64}
   error_variances: Judgement{Int64}
   error_covariances_within: Judgement{Int64}
   error_covariances_between: Judgement{Int64}
   crossloadings_incoming: Judgement{Int64}
   crossloadings_outgoing: Judgement{Int64}
```
"""
#composite type
struct CFAFactor <: AbstractCFA
    type1::SimpleCFA
    type2::HierarchicalCFA
end

struct SimpleCFA <: AbstractCFA
    n_sample::Judgement{ <: Union{ <:Int, Missing}}
    n_variables::Judgement{ <: Union{ <:Int, Missing}}
    loadings::Judgement{ <: Union{ <: AbstractArray{ <: Number}, Missing}}
    factor_variance::Judgement{ <: Union{ <:Number, Missing}}
    error_variances::Judgement{<:Union{<:AbstractArray{<:Number},<: Int, Missing}}
    error_covariances_within::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    error_covariances_between::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    crossloadings_incoming::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    crossloadings_outgoing::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    SimpleCFA(n_sample, n_variables, loadings,factor_variance, error_variances, error_covariances_within, error_covariances_between, crossloadings_incoming, crossloadings_outgoing) =
        new(J(n_sample), J(n_variables), J(loadings), J(factor_variance), J(error_variances), J(error_covariances_within), J(error_covariances_between), J(crossloadings_incoming), J(crossloadings_outgoing))
end

# this struct needs fields to specify for second order factors. And hanlde latent covariances
struct HierarchicalCFA <: AbstractCFA
    n_sample::Judgement{ <: Union{ <:Int, Missing}}
    n_variables::Judgement{ <: Union{ <:Int, Missing}}
    loadings::Judgement{ <: Union{ <: AbstractArray{ <: Number}, Missing}}
    factor_variance::Judgement{ <: Union{ <:Number, Missing}}
    error_variances::Judgement{<:Union{<:AbstractArray{<:Number},<: Int, Missing}}
    error_covariances_within::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    error_covariances_between::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    crossloadings_incoming::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    crossloadings_outgoing::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    HierarchicalCFA(n_sample, n_variables, loadings,factor_variance, error_variances, error_covariances_within, error_covariances_between, crossloadings_incoming, crossloadings_outgoing) =
        new(J(n_sample), J(n_variables), J(loadings), J(factor_variance), J(error_variances), J(error_covariances_within), J(error_covariances_between), J(crossloadings_incoming), J(crossloadings_outgoing))
end

#multiple dispatch function for CFAFactor; insert either simple or hierachical CFA
function CFAFactor(;n_sample = missing,
        n_variables,
        loadings,
        factor_variance,
        error_variances = 0,
        error_covariances_within = 0,
        error_covariances_between = 0,
        crossloadings_incoming = 0,
        crossloadings_outgoing = 0)

    simple_cfa = SimpleCFA(n_sample, n_variables, loadings, factor_variance, error_variances, 
           error_covariances_within, error_covariances_between, crossloadings_incoming, 
           crossloadings_outgoing)
    # here also adjust for changes to come in struct above.
    hierarchical_cfa = HierarchicalCFA(n_sample, n_variables, loadings, factor_variance, error_variances, 
                       error_covariances_within, error_covariances_between, crossloadings_incoming, 
                       crossloadings_outgoing)

return CFAFactor(simple_cfa, hierarchical_cfa)
end



"""
AbstractCFA Taxons.
Consists of Factors (measurement model) and a graph from StenoGraphs (structural model). 

## Arguments

- `measurement_model`: Vector of Factors.
- `structural_model`: Graph from StenoGraphs package. Defines the latent relations between the factors of measurement_model.  

```julia
using StenoGraphs
using Taxonomy

factor1 = CFAFactor(n_variables = 2, loadings = [1, 0.4], factor_variance = 0.7)
factor2 = CFAFactor(n_variables = 2, loadings = [0.7, 0.3], factor_variance = 1)

graph = @StenoGraph begin
    # latent regressions
    fac1 → fac2
end

AbstractCFA(measurement_model = [factor1, factor2], 
structural_model = graph)

# output

AbstractCFA
   n_sample: Judgement{Missing}
   measurement_model: Judgement{Vector{Factor}}
   structural_model: Judgement{Vector{DirectedEdge{SimpleNode{Symbol}, SimpleNode{Symbol}}}}

```
"""
struct HierarchicalCFA <: Taxon
    n_sample::Judgement{ <: Union{ <:Int, Missing}}
    measurement_model::Judgement{ <: Union{<:AbstractArray{<: Factor}, Missing}}
    structural_model::Judgement{ <: Union{<:AbstractArray{<: StenoGraphs.AbstractEdge}, Missing}}
    CFA(n_sample, measurement_model, structural_model) = 
    new(J(n_sample), J(measurement_model), J(structural_model))
end

struct SimpleCFA <: Taxon
    n_sample::Judgement{ <: Union{ <:Int, Missing}}
    measurement_model::Judgement{ <: Union{<:AbstractArray{<: Factor}, Missing}}
    hierarchical_model::Judgement{ <: Union{<:AbstractArray{<: StenoGraphs.AbstractEdge}, Missing}}
    CFA(n_sample, measurement_model, structural_model) = 
    new(J(n_sample), J(measurement_model), J(structural_model))
end

function CFA(;n_sample = missing,
    measurement_model,
    hierarchical_model) # hier hierarchical CFA model festlegen
    CFA(n_sample, measurement_model, hierarchical_model)
end

