abstract type AbstractLocation end
abstract type AbstractDOI <: AbstractLocation end
abstract type AbstractMeta end
"""
Taxon is the supertype of all taxons.
"""

abstract type Taxon end
abstract type NoTaxon end
abstract type CFA <: Taxon end
abstract type AbstractSEM <: Taxon end


struct Pathmodel <: Taxon end
struct SimpleCFA <: CFA end
struct HierarchicalCFA <: CFA end
struct CrossSectionalSEM <: SEM end
struct LongitudinalSEM <: SEM end

#composite:
struct CFAFactor <: CFA end
    type1::SimpleCFA
    type2::HierarchicalCFA

struct SEMFactor <: AbstractSEM end
    type1::CrossSectionalSEM
    type2::LongitudinalSEM

# multiple dispatch:
function CFAFactor(factor::CFAFactor)
    # do something with SimpleCFA and HierarchicalCFA
    println("Doing something with CFAFactor")
end

function SEMFactor(factor::SEMSector)
    # do something with CrossSectionalSEM and LongitudinalSEM
    println("Doing something with SEMFactor")
end
