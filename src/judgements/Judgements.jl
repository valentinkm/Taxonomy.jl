module Judgements
    import Taxonomy: AbstractJudgement, Taxon, JudgementLevel
    check_certainty(c) = ((c < 0.0) || (c > 1.0)) ? throw(ArgumentError("Certainty must be between 0 and 1.")) : nothing
    export AnyLevelJudgement, RecordJudgement, StudyJudgement, ModelJudgement
    export check_judgement_level, correct_judgement_level, judgements
    include("level.jl")

    include("fixed.jl")

    export @newjudgement
    export J, Judgement, NoJudgement, convert, rating, comment, certainty, judgement_key
    include("judgement.jl")

    # exports via code gen within the file
    include("constrained.jl")
    
    export judgement_dict
    include("dict.jl")

    export CFI, Empirical, Lang, N
    include("predefined_judgements.jl")
end
