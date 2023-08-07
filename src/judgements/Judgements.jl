module Judgements
    import Taxonomy: AbstractJudgement
    check_certainty(c) = ((c < 0.0) || (c > 1.0)) ? throw(ArgumentError("Certainty must be between 0 and 1.")) : nothing
    export J, Judgement, NoJudgement, convert, rating, location, certainty
    include("judgement.jl")
    export JudgementBool, JudgementNumber, JudgementInt, JudgementFloat, JudgementString
    export @newjudgement
    include("constrained.jl")
    export JudgementLanguage
    include("language.jl")
end