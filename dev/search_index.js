var documenterSearchIndex = {"docs":
[{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"It is often the case that you can not code every information we would like to have. Principly, you are not required to write anything:","category":"page"},{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"DocTestFilters = [r\"May we suggest: .*\", r\"└ @ Taxonomy .*\"]","category":"page"},{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"julia> using Taxonomy\n\njulia> Record()\n┌ Warning: Please provide your rater ID. This should be your initials.\n└ @ Taxonomy ~/work/Taxonomy.jl/Taxonomy.jl/src/record.jl:15\n┌ Warning: You really should supply an ID. May we suggest: 5260937a-209d-4cf9-a5d5-120f350fae45\n└ @ Taxonomy ~/work/Taxonomy.jl/Taxonomy.jl/src/record.jl:14\n┌ Warning: You really should supply a location.\n└ @ Taxonomy ~/work/Taxonomy.jl/Taxonomy.jl/src/record.jl:19\n┌ Warning: Some of the metadata seem to be incomplete. Check again.\n└ @ Taxonomy ~/work/Taxonomy.jl/Taxonomy.jl/src/record.jl:25\n┌ Warning: `taxons` is missing. Maybe you mean `NoTaxon()`?\n└ @ Taxonomy ~/work/Taxonomy.jl/Taxonomy.jl/src/record.jl:29\n┌ Warning: `spec` is missing. Maybe you mean `NoJudgment()`?\n└ @ Taxonomy ~/work/Taxonomy.jl/Taxonomy.jl/src/record.jl:32\n┌ Warning: `data` is missing. Maybe you mean `NoJudgment()`?\n└ @ Taxonomy ~/work/Taxonomy.jl/Taxonomy.jl/src/record.jl:35\nRecord(missing, missing, NoLocation(), IncompleteMeta(missing, missing, missing), missing, missing, missing)","category":"page"},{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"DocTestFilters = nothing","category":"page"},{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"As you probably notice we warn you to do that. This is to encourage you to think twice, however, after having thought twice about it, you may silence every warning with explicitly suppliyng \"empty\" instances (except ID, really nothing should hinder you to supply a random id).","category":"page"},{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"Record(rater = \"AP\",\nid = \"6ca721fe-619e-42cc-ad8b-047c5e0451e5\",\nlocation = NoLocation(),\nmeta = MetaData(missing, missing, missing),\ntaxons = [NoTaxon()],\nspec = NoJudgement(),\ndata = NoJudgement())\n\n# output\n\nRecord(\"AP\", UUID(\"6ca721fe-619e-42cc-ad8b-047c5e0451e5\"), NoLocation(), IncompleteMeta(missing, missing, missing), [NoTaxon()], Judgement{Missing}(missing, 0.0, missing), Judgement{Missing}(missing, 0.0, missing))","category":"page"},{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"We therefore differenciate between \"lazy\" missings and intentional missings. The former remind you that you missed them, the latter will not bother you.","category":"page"},{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"Another, unfortunately more complicated distinction, is between things you have tried to find out and determined that they probably are not accessable, and things that you havened made your mind up about. We assume that missings in about metadata belong to the first category (you tried to find out, but could not) and that missings about judgements about the paper belong to the second category (you have not bothered to find out). We make this distinction because you might think that it takes to much time to find something out and want to express that you are not certain it exists. But sometimes you have checked but it really does not exist. E.g. at first glance you have not found a dataset, than use NoJudgement(). It implies that you are absolutly uncertain that this is missing:","category":"page"},{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"julia> NoJudgement()\nJudgement{Missing}(missing, 0.0, missing)","category":"page"},{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"Or you have checked everywhere but there does not seem to be any data, than use:","category":"page"},{"location":"tutorials/missing/","page":"Something is missing","title":"Something is missing","text":"julia> Judgement(false, 1.0) # false = no data, 1 = certain\nJudgement{Bool}(false, 1.0, missing)","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Taxonomy","category":"page"},{"location":"#Taxonomy","page":"Home","title":"Taxonomy","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Taxonomy.","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Taxons","page":"Home","title":"Taxons","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Taxon\nFactor\nCFA","category":"page"},{"location":"#Taxonomy.Taxon","page":"Home","title":"Taxonomy.Taxon","text":"Taxon is the supertype of all taxons.\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.Factor","page":"Home","title":"Taxonomy.Factor","text":"Factor Taxon.  Building Block for CFA Taxonomy. Multiple Factors can be combined to a CFA. \n\nArguments\n\nn_sample: Number of observed cases. May be between different taxons from same paper sometime, e.g. multigroup models.\nn_variables: Number of variables (possibly observed/manifest).\nloadings: Vector of loadings, one for each item. \nerror_covariances_within: Vector of covariances within factor.\nerror_covariances_between: Vector of covariances the factor shares with a different factor. \ncrossloadings_incoming: Vector of crossloadings coming from other factors. They should be lower than the loading coming to the item from this factor.  \ncrossloadings_outgoing: Vector of crossloadings going to other items which have higher loadings from other factors. \n\nFactor(n_variables = 2, loadings = [1, 0.4])\n\n# output\nFactor(Judgement{Missing}(missing, 1.0, missing), Judgement{Int64}(2, 1.0, missing), Judgement{Vector{Float64}}([1.0, 0.4], 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing))\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.CFA","page":"Home","title":"Taxonomy.CFA","text":"CFA Taxon.  Consists of Factors (measurement model) and a graph from StenoGraphs (structural model). \n\nArguments\n\nmeasurement_model: Vector of Factors.\nstructural_model: Graph from StenoGraphs package. Defines the latent relations between the factors of measurement_model.  \n\nusing StenoGraphs\nusing Taxonomy\n\nfactor1 = Factor(n_variables = 2, loadings = [1, 0.4])\nfactor2 = Factor(n_variables = 2, loadings = [0.7, 0.3])\n\ngraph = @StenoGraph begin\n    # latent regressions\n    fac1 → fac2\nend\n\nCFA(measurement_model = [factor1, factor2], \nstructural_model = graph )\n\n# output\nCFA(Judgement{Missing}(missing, 1.0, missing), Judgement{Vector{Factor}}(Factor[Factor(Judgement{Missing}(missing, 1.0, missing), Judgement{Int64}(2, 1.0, missing), Judgement{Vector{Float64}}([1.0, 0.4], 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing)), Factor(Judgement{Missing}(missing, 1.0, missing), Judgement{Int64}(2, 1.0, missing), Judgement{Vector{Float64}}([0.7, 0.3], 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing))], 1.0, missing), Judgement{Vector{DirectedEdge{SimpleNode{Symbol}, SimpleNode{Symbol}}}}(fac1 → fac2\n, 1.0, missing))\n\n\n\n\n\n","category":"type"},{"location":"#Extractors","page":"Home","title":"Extractors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"n_sample","category":"page"},{"location":"#Taxonomy.n_sample","page":"Home","title":"Taxonomy.n_sample","text":"Function to extract sample size.\n\nArguments\n\nx: Something of type Taxon\n\nReturn\n\nReturns a Judgement\n\nf = Factor(n_sample = 100, n_variables = 2, loadings = [1, 0.4])\nrating(n_sample(f))\n\n# output\n100\n\n\n\n\n\n","category":"function"},{"location":"#ID","page":"Home","title":"ID","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"generate_id","category":"page"},{"location":"#Taxonomy.generate_id","page":"Home","title":"Taxonomy.generate_id","text":"Generate an Entry ID\n\nTo create links between entries we need a stable reference point. This ID is generated initially from url(location) and if the url is missing, it is generated randomly. After the ID is generated once, it is saved with the Record and should not be changed.\n\n\n\n\n\n","category":"function"},{"location":"#Judgement","page":"Home","title":"Judgement","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Judgement\nJ\nNoJudgement\nrating\nlocation(::Judgement)\ncertainty","category":"page"},{"location":"#Taxonomy.Judgement","page":"Home","title":"Taxonomy.Judgement","text":"A judgment about any parameter etc.\n\nArguments\n\nrating: The rating, e.g. \"CFA\" or 1.0.\ncertainty: If uncertain, a number between 0.0 and 1.0 (0-100%)\nlocation: optional, Where in the Paper PDF was the location retieved, e.g. section, page, table number, figure number.\n\njulia> Judgement(1.0, .99, \"Figure 1\");\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.J","page":"Home","title":"Taxonomy.J","text":"Shorthand for Judgement\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.NoJudgement","page":"Home","title":"Taxonomy.NoJudgement","text":"Abstaining from any judgement.\n\nThis implies that your best guess is missing and you are absolutely uncertain about this judgement.\n\njulia> NoJudgement()\nJudgement{Missing}(missing, 0.0, missing)\n\n\n\n\n\n","category":"function"},{"location":"#Taxonomy.rating","page":"Home","title":"Taxonomy.rating","text":"Extract rating from Judgement.\n\nIf rating is called on a Judgement it returns the rating, on everything it returns identity.\n\n\n\n\n\n","category":"function"},{"location":"#Taxonomy.location-Tuple{Judgement}","page":"Home","title":"Taxonomy.location","text":"Extract location from Judgement.\n\n\n\n\n\n","category":"method"},{"location":"#Taxonomy.certainty","page":"Home","title":"Taxonomy.certainty","text":"Extract certainty from Judgement.\n\n\n\n\n\n","category":"function"},{"location":"#Metadata","page":"Home","title":"Metadata","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MetaData\napa\njson\nauthor\nyear\njournal\nMinimalMeta\nIncompleteMeta\nExtensiveMeta","category":"page"},{"location":"#Taxonomy.MetaData","page":"Home","title":"Taxonomy.MetaData","text":"Save metadata.\n\nCan be from complete minimal metadata, incomplete metadata or preferably from DOI.\n\njulia> min = MetaData(\"Peikert, Aaron\", 2022, \"Journal of Statistical Software\");\n\njulia> incomplete = MetaData(\"Peikert, Aaron\", 2022, missing);\n\njulia> extensive = MetaData(DOI(\"10.5281/zenodo.6719627\"));\n\n\n\n\n\n\n","category":"function"},{"location":"#Taxonomy.apa","page":"Home","title":"Taxonomy.apa","text":"Get an APA citation.\n\njulia> apa(DOI(\"10.5281/zenodo.6719627\"))\n\"Ernst, M. S., &amp; Peikert, A. (2022). <i>StructuralEquationModels.jl</i> (Version v0.1.0) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.6719627\"\n\n\n\n\n\n","category":"function"},{"location":"#Taxonomy.json","page":"Home","title":"Taxonomy.json","text":"Get a Citeproc JSON.\n\nCSL JSON Documentation\n\nCSL JSON can be read by Zotero and automatically generated by doi.org from DOI. All availible information are included and saved in a Dict.\n\njulia> json(DOI(\"10.5281/zenodo.6719627\"))\nDict{String, Any} with 11 entries:\n  \"publisher\" => \"Zenodo\"\n  \"issued\"    => Dict{String, Any}(\"date-parts\"=>Any[Any[2022, 6, 24]])\n  \"author\"    => Any[Dict{String, Any}(\"family\"=>\"Ernst\", \"given\"=>\"Maximilian …\n  \"id\"        => \"https://doi.org/10.5281/zenodo.6719627\"\n  \"copyright\" => \"MIT License\"\n  \"version\"   => \"v0.1.0\"\n  \"DOI\"       => \"10.5281/ZENODO.6719627\"\n  \"URL\"       => \"https://zenodo.org/record/6719627\"\n  \"title\"     => \"StructuralEquationModels.jl\"\n  \"abstract\"  => \"StructuralEquationModels v0.1.0 This is a package for Structu…\n  \"type\"      => \"book\"\n\n\n\n\n\n","category":"function"},{"location":"#Taxonomy.author","page":"Home","title":"Taxonomy.author","text":"Extract the author.\n\njulia> doi = MetaData(DOI(\"10.1126/SCIENCE.169.3946.635\"));\n\njulia> author(doi)\n\"Frank, Henry S.\"\n\n\n\n\n\n","category":"function"},{"location":"#Taxonomy.year","page":"Home","title":"Taxonomy.year","text":"Extract the year.\n\njulia> doi = MetaData(DOI(\"10.1126/SCIENCE.169.3946.635\"));\n\njulia> year(doi)\n1970\n\n\n\n\n\n","category":"function"},{"location":"#Taxonomy.journal","page":"Home","title":"Taxonomy.journal","text":"Extract the journal.\n\njulia> doi = MetaData(DOI(\"10.1126/SCIENCE.169.3946.635\"));\n\njulia> journal(doi)\n\"Science\"\n\n\n\n\n\n","category":"function"},{"location":"#Taxonomy.MinimalMeta","page":"Home","title":"Taxonomy.MinimalMeta","text":"A representation of the most important metadata.\n\njulia> min = MetaData(\"Peikert, Aaron\", 2022, \"Journal of Statistical Software\");\n\njulia> typeof(min)\nMinimalMeta\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.IncompleteMeta","page":"Home","title":"Taxonomy.IncompleteMeta","text":"A representation of Metadata when we can not even capture the most important metadata.\n\njulia> incomplete = MetaData(missing, 2022, \"Journal of Statistical Software\");\n\njulia> typeof(incomplete)\nIncompleteMeta\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.ExtensiveMeta","page":"Home","title":"Taxonomy.ExtensiveMeta","text":"The metadata we can gather from doi.org.\n\njulia> doi = MetaData(DOI(\"10.1126/SCIENCE.169.3946.635\"));\n\njulia> typeof(doi)\nExtensiveMeta{MinimalMeta}\n\n\n\n\n\n","category":"type"},{"location":"#DOI","page":"Home","title":"DOI","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DOI\nUsualDOI\nUnusualDOI\nNoDOI\nNoLocation\nurl\nvalid_doi","category":"page"},{"location":"#Taxonomy.DOI","page":"Home","title":"Taxonomy.DOI","text":"Alias for UsualDOI.\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.UsualDOI","page":"Home","title":"Taxonomy.UsualDOI","text":"Construct a validated DOI\n\nMost valid DOIs (not all) can be simply validated via a regular expression.\n\nArguments\n\ndoi::String: a DOI without resolver (e.g. without doi.org), capitalization does not matter\nfallback::String: an optional fallback link where one maybe can find the content in case the doi fails\n\njulia> DOI(\"10.5281/zenodo.6719627\")\nUsualDOI{String, Missing}(\"10.5281/ZENODO.6719627\", missing)\n\njulia> DOI(\"10.5281/zenodo.6719627\", \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\")\nUsualDOI{String, String}(\"10.5281/ZENODO.6719627\", \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\")\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.UnusualDOI","page":"Home","title":"Taxonomy.UnusualDOI","text":"Construct an unvalidated DOI\n\nYou should prefer an validated UsualDOI but if you have tested the DOI and are sure it links were it supposed to link, go ahead and create an unvalidated doi.\n\njulia> UnusualDOI(\"weird10.5281doi/zenodo.6719627\")\nUnusualDOI{String, Missing}(\"WEIRD10.5281DOI/ZENODO.6719627\", missing)\n\njulia> UnusualDOI(\"weird10.5281doi/zenodo.6719627\", \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\")\nUnusualDOI{String, String}(\"WEIRD10.5281DOI/ZENODO.6719627\", \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\")\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.NoDOI","page":"Home","title":"Taxonomy.NoDOI","text":"What to do if there is no doi\n\nLast resort if there is no DOI. Than we save other metadata, similar to BibTex.\n\nArguments\n\nurl::String: an link where one maybe can find the content in case the doi fails\nauthor::String: like in BibTex, e.g. \"Peikert, Aaron and Ernst, Maximilian S. and Bode, Clifford\"\ndate::Union{Date, Missing}: optional date\nyear::Union{Int64}: optional if date is supplied\njournal::String: The outlet of the publication\nother::Dict: more BibTexlike metadata\n\nNoDOI(\n    url = \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\",\n    author = \"Ernst, Maximilian Stefan and Peikert, Aaron\",\n    title = \"StructuralEquationModels.jl: A fast and flexible SEM framework\",\n    date = Date(\"2022-06-24\"), # year is inferred\n    journal = \"No Real Journal\",\n    awesome = \"Yes\", # other metadata\n    software = \"naturally\", # some more metadata\n    citations = 500\n)\nNoDOI(\n    url = \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\",\n    author = \"Ernst, Maximilian Stefan and Peikert, Aaron\",\n    title = \"StructuralEquationModels.jl: A fast and flexible SEM framework\",\n    year = 2022, # date is omitted\n    journal = \"No Real Journal\"\n)\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.NoLocation","page":"Home","title":"Taxonomy.NoLocation","text":"When everything fails.\n\nThis is a placeholder if really no location can be found.\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.url","page":"Home","title":"Taxonomy.url","text":"Get URL from location.\n\njulia> url(DOI(\"10.1126/SCIENCE.169.3946.635\"))\n\"https://doi.org/10.1126/SCIENCE.169.3946.635\"\n\nlocation = NoDOI(\n    url = \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\",\n    author = \"Ernst, Maximilian Stefan and Peikert, Aaron\",\n    title = \"StructuralEquationModels.jl: A fast and flexible SEM framework\",\n    year = 2022, # date is omitted\n    journal = \"No Real Journal\"\n)\n\nurl(location)\n\n# output\n\n\"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\"\n\n\n\n\n\n","category":"function"},{"location":"#Taxonomy.valid_doi","page":"Home","title":"Taxonomy.valid_doi","text":"Validate DOI via Regex\n\nRegular expression taken from:\n\nhttps://www.crossref.org/blog/dois-and-matching-regular-expressions/\n\n\n\n\n\n","category":"function"},{"location":"tutorials/start/","page":"Getting Started","title":"Getting Started","text":"The idea behind doing the coding in Julia is that we can automate two things.","category":"page"},{"location":"tutorials/start/","page":"Getting Started","title":"Getting Started","text":"We can enforce syntactical requirements with a quick feedback loop for the coder, resulting in less errors.\nWe can automatically infer much information, resulting in speedier and less error prone coding.","category":"page"},{"location":"tutorials/start/","page":"Getting Started","title":"Getting Started","text":"We implement these automation to detect error/misunderstanding/etc while coding and not long after.","category":"page"},{"location":"tutorials/start/","page":"Getting Started","title":"Getting Started","text":"To get a feeling for how this is supposed to work, we code the famous political democracy paper together:","category":"page"},{"location":"tutorials/start/","page":"Getting Started","title":"Getting Started","text":"DocTestFilters = [r\"\\\".*\\\"\", r\"May we suggest: .*\", r\"└ @ Taxonomy .*\"]","category":"page"},{"location":"tutorials/start/","page":"Getting Started","title":"Getting Started","text":"using Taxonomy \nfirst_record = Record(\n    rater = \"AP\",\n    location = DOI(\"10.2307/2095172\"),\n    taxons = [Factor(n_variables = 6, loadings = [1, 1.19, 0.53, 0.91, 1, 1], error_covariances_within = [10.7, 12.9, 19])],\n    spec = true,\n    data = true\n)\n\n# output\n\n┌ Warning: You really should supply an ID. May we suggest (from DOI): 8f1713c9-482b-58cb-8ed4-128c03e9dafb\n└ @ Taxonomy ~/Downloads/Taxonomy.jl/src/record.jl:19\nRecord(\"AP\", missing, UsualDOI{String, Missing}(\"10.2307/2095172\", missing), ExtensiveMeta{MinimalMeta}(MinimalMeta(\"Bollen, Kenneth A.\", 1980, \"American Sociological Review\"), \"Bollen, K. A. (1980). Issues in the Comparative Measurement of Political Democracy. American Sociological Review, 45(3), 370. https://doi.org/10.2307/2095172\\n\", Dict{String, Any}(\"publisher\" => \"SAGE Publications\", \"member\" => \"179\", \"issue\" => \"3\", \"indexed\" => Dict{String, Any}(\"date-parts\" => Any[Any[2023, 1, 12]], \"date-time\" => \"2023-01-12T14:50:52Z\", \"timestamp\" => 1673535052942), \"reference-count\" => 0, \"container-title-short\" => \"American Sociological Review\", \"issued\" => Dict{String, Any}(\"date-parts\" => Any[Any[1980, 6]]), \"container-title\" => \"American Sociological Review\", \"subject\" => Any[\"Sociology and Political Science\"], \"resource\" => Dict{String, Any}(\"primary\" => Dict{String, Any}(\"URL\" => \"http://www.jstor.org/stable/2095172?origin=crossref\"))…)), Factor[Factor(Judgement{Missing}(missing, 1.0, missing), Judgement{Int64}(6, 1.0, missing), Judgement{Vector{Float64}}([1.0, 1.19, 0.53, 0.91, 1.0, 1.0], 1.0, missing), Judgement{Vector{Float64}}([10.7, 12.9, 19.0], 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing))], Judgement{Bool}(true, 1.0, missing), Judgement{Bool}(true, 1.0, missing))\n","category":"page"},{"location":"tutorials/start/","page":"Getting Started","title":"Getting Started","text":"DocTestFilters = nothing","category":"page"},{"location":"tutorials/start/","page":"Getting Started","title":"Getting Started","text":"julia> year(first_record)\n1980\njulia> apa(first_record)\n\"Bollen, K. A. (1980). Issues in the Comparative Measurement of Political Democracy. American Sociological Review, 45(3), 370. https://doi.org/10.2307/2095172\\n\"\njulia> spec(first_record)\nJudgement{Bool}(true, 1.0, missing)\njulia> taxons(first_record)\n1-element Vector{Factor}:\n Factor(Judgement{Missing}(missing, 1.0, missing), Judgement{Int64}(6, 1.0, missing), Judgement{Vector{Float64}}([1.0, 1.19, 0.53, 0.91, 1.0, 1.0], 1.0, missing), Judgement{Vector{Float64}}([10.7, 12.9, 19.0], 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing), Judgement{Int64}(0, 1.0, missing))","category":"page"},{"location":"tutorials/start/","page":"Getting Started","title":"Getting Started","text":"Generally, DOIs are the best and easiest thing to get metadata, however, sometimes none is availible:","category":"page"},{"location":"tutorials/start/","page":"Getting Started","title":"Getting Started","text":"Record(\n    location = NoDOI(\n        url = \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\",\n        author = \"Ernst, Maximilian Stefan and Peikert, Aaron\",\n        title = \"StructuralEquationModels.jl: A fast and flexible SEM framework\",\n        year = 2022, # date is omitted\n        journal = \"No Real Journal\"\n    ),\n    taxons = [Factor(n_variables = 6, loadings = [1, 1.19, 0.53, 0.91, 1, 1], error_covariances_within = [10.7, 12.9, 19])],\n    spec = true,\n    data = true\n)","category":"page"}]
}
