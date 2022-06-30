var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Taxonomy","category":"page"},{"location":"#Taxonomy","page":"Home","title":"Taxonomy","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Taxonomy.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Taxonomy]","category":"page"},{"location":"#Taxonomy.DOI","page":"Home","title":"Taxonomy.DOI","text":"Alias for UsualDOI.\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.UnusualDOI","page":"Home","title":"Taxonomy.UnusualDOI","text":"Construct an unvalidated DOI\n\nYou should prefer an validated UsualDOI but if you have tested the DOI and are sure it links were it supposed to link, go ahead and create an unvalidated doi.\n\njulia> UnusualDOI(\"weird10.5281doi/zenodo.6719627\")\nUnusualDOI{String, Missing}(\"WEIRD10.5281DOI/ZENODO.6719627\", missing)\n\njulia> UnusualDOI(\"weird10.5281doi/zenodo.6719627\", \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\")\nUnusualDOI{String, String}(\"WEIRD10.5281DOI/ZENODO.6719627\", \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\")\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.UsualDOI","page":"Home","title":"Taxonomy.UsualDOI","text":"Construct a validated DOI\n\nMost valid DOIs (not all) can be simply validated via a regular expression.\n\nArguments\n\nobserved::String: a DOI without resolver (e.g. without doi.org), capitalization does not matter\nfallback::String: an optional fallback link where one maybe can find the content in case the doi fails\n\njulia> DOI(\"10.5281/zenodo.6719627\")\nUsualDOI{String, Missing}(\"10.5281/ZENODO.6719627\", missing)\n\njulia> DOI(\"10.5281/zenodo.6719627\", \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\")\nUsualDOI{String, String}(\"10.5281/ZENODO.6719627\", \"https://github.com/StructuralEquationModels/StructuralEquationModels.jl\")\n\n\n\n\n\n","category":"type"},{"location":"#Taxonomy.valid_doi-Tuple{Any}","page":"Home","title":"Taxonomy.valid_doi","text":"Validate DOI via Regex\n\nRegular expression taken from:\n\nhttps://www.crossref.org/blog/dois-and-matching-regular-expressions/\n\n\n\n\n\n","category":"method"}]
}
