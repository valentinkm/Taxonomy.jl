abstract type AbstractLocation end
abstract type AbstractDOI <: AbstractLocation end
abstract type AbstractMeta end
"""
Taxon is the supertype of all taxons.
"""

abstract type Taxon end
abstract type NoTaxon end
abstract type AbstractCFA <: Taxon end
abstract type AbstractSEM <: Taxon end


struct Pathmodel <: Taxon end
struct SimpleCFA <: AbstractCFA end
struct HierarchicalCFA <: AbstractCFA end
struct CrossSectionalSEM <: AbstractSEM end
struct LongitudinalSEM <: AbstractSEM end

#composite:
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
