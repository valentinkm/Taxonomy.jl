abstract type AbstractLocation end
abstract type AbstractDOI <: AbstractLocation end
abstract type AbstractMeta end

"""
Factor is the supertype for all measurement models:
"""
abstract type AbstractFactor end
# this does not belong here, needs to be moved to something like factor.jl
# this struct needs all possible fields of all possible factors. e.g. time coding LGCM
struct Factor <: AbstractFactor
    n_sample::Judgement{ <: Union{ <:Int, Missing}}
    n_variables::Judgement{ <: Union{ <:Int, Missing}}
    loadings::Judgement{ <: Union{ <: AbstractArray{ <: Number}, Missing}}
    factor_variance::Judgement{ <: Union{ <:Number, Missing}}
    error_variances::Judgement{<:Union{<:AbstractArray{<:Number},<: Int, Missing}}
    error_covariances_within::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    error_covariances_between::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    crossloadings_incoming::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    crossloadings_outgoing::Judgement{ <: Union{ <: AbstractArray{ <: Number}, <: Int, Missing}}
    Factor(n_sample, n_variables, loadings,factor_variance, error_variances, error_covariances_within, error_covariances_between, crossloadings_incoming, crossloadings_outgoing) =
        new(J(n_sample), J(n_variables), J(loadings), J(factor_variance), J(error_variances), J(error_covariances_within), J(error_covariances_between), J(crossloadings_incoming), J(crossloadings_outgoing))
end

struct PathFactor <: AbstractFactor
    n_sample::Judgement{ <: Union{ <:Int, Missing}}
    factor_variance::Judgement{ <: Union{ <:Number, Missing}}
    PathFactor(n_sample, factor_variance) =
        new(J(n_sample), J(factor_variance))
end



"""
Taxon is the supertype of all taxons.
"""
abstract type Taxon end
abstract type NoTaxon end

abstract type AbstractCFA <: Taxon end
abstract type AbstractSEM <: Taxon end


struct Pathmodel <: Taxon end

struct CrossSectionalSEM <: AbstractSEM end
struct LongitudinalSEM <: AbstractSEM end

"""
#composite:
"""
struct TaxonFactor
    taxon::Taxon
    factor::AbstractFactor
end


"""
struct CFAFactor <: AbstractCFA
    type1::SimpleCFA
    type2::HierarchicalCFA
end

struct SEMFactor <: AbstractSEM
    type1::CrossSectionalSEM
    type2::LongitudinalSEM
end

# multiple dispatch:
function CFAFactor(factor::CFAFactor)
    # do something with SimpleCFA and HierarchicalCFA
    println("Doing something with CFAFactor")
end

function SEMFactor(factor::SEMFactor)
    # do something with CrossSectionalSEM and LongitudinalSEM
    println("Doing something with SEMFactor")
end

"""