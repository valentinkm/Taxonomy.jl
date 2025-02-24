function invariant_pair(x::Pair{K,V}) where {K,V<:Record}
    x[1] == id(x[2]) ? x : throw(ArgumentError("Key and Record id must be the same."))
end

function invariant_pair(key, value::Record)
    invariant_pair(Pair(key, value))
end

function invariant_pair(value::Record)
    invariant_pair(id(value), value)
end


"""
RecordDatabase(records::Record...)

`Record`s need to be stored somewhere.

```jldoctest
julias> RecordDatabase()
RecordDatabase{Base.UUID, Record}()

julia> first = Record(rater = "AP", id = "552ef675-5c7b-4ce1-880b-c45b833fdfcb", location = NoLocation(), meta = MetaData(missing, missing, missing), judgements = Dict("Taxon" => J(NoTaxon())), spec = NoJudgement(), data = NoJudgement());

julias> second = Record(rater = "AP", id = "58c55701-0362-40c7-849c-5d12e5026238", location = NoLocation(), meta = MetaData(missing, missing, missing), judgements = Dict("Taxon" => J(NoTaxon())), spec = NoJudgement(), data = NoJudgement());

julia> rd = RecordDatabase(first, second)
RecordDatabase{Base.UUID, Record} with 2 entries:
  UUID("58c55701-0362-40c7-849c-5d12e5026238") => Record…
  UUID("552ef675-5c7b-4ce1-880b-c45b833fdfcb") => Record…

julias> rd += Record(rater = "AP", id = "2d7ad584-8ec0-47dd-807b-280fba2978f8", location = NoLocation(), meta = MetaData(missing, missing, missing), judgements = Dict("Taxon" => J(NoTaxon())), spec = NoJudgement(), data = NoJudgement())
RecordDatabase{Base.UUID, Record} with 3 entries:
  UUID("2d7ad584-8ec0-47dd-807b-280fba2978f8") => Record…
  UUID("58c55701-0362-40c7-849c-5d12e5026238") => Record…
  UUID("552ef675-5c7b-4ce1-880b-c45b833fdfcb") => Record…
```
"""

struct RecordDatabase{K<:UUID,V<:Record} <: Base.AbstractDict{K,V}
    records::Dict{K,V}
end

@inline Base.length(rd::RecordDatabase) = length(rd.records)
@inline Base.iterate(rd::RecordDatabase) = Base.iterate(rd.records)
@inline Base.iterate(rd::RecordDatabase, state) = Base.iterate(rd.records, state)


# Constructor to create an empty RecordDatabase
RecordDatabase() = RecordDatabase(Dict{UUID,Record}())
RecordDatabase(records::Record...) = RecordDatabase(Dict(invariant_pair.(records)))
Base.convert(::Type{Pair{UUID,Record}}, r::Record) = invariant_pair(r)

# Overload the push! function to add a new Record to the RecordDatabase
# builds function to be able to push new records to the database
# push! as function with two arguments: db of type RecordDatabase and record of type Record, returns an instance of RecordDatabase
function Base.push!(x::RecordDatabase, new::Pair)
    check_uuid(x, new)
    check_url(x, new)
    push!(x.records, invariant_pair(new))
    x
end
Base.push!(x::RecordDatabase{K,V}, new) where {K,V} = push!(x, convert(Pair{K,V}, new))
Base.setindex!(rd::RecordDatabase{K,V}, value::V, key::K) where {K,V} =
    push!(rd, invariant_pair(key, value))
Base.haskey(rd::RecordDatabase{K,V}, key::K) where {K,V} = haskey(rd.records, key)
Base.get(rd::RecordDatabase, args...) = get(rd.records, args...)
function Base.delete!(rd::RecordDatabase{K,V}, key::K) where {K,V}
    delete!(rd.records, key)
    rd
end
function Base.merge(x::RecordDatabase, y::RecordDatabase)
    check_uuid(x, y)
    check_url(x, y)
    return RecordDatabase(merge(x.records, y.records))
end
Base.:+(x::RecordDatabase, y::RecordDatabase) = merge(x, y)
Base.:+(x::RecordDatabase, y::Record) = push!(x, y)
Base.:+(x::Record, y::RecordDatabase) = y + x


# Base.setindex!(rd::RecordDatabase{K, V}, value::K, key::V) where {K, V} = push!(rd, key => value)
# Base.haskey(rd::RecordDatabase{K, V}, key::K) where {K, V} = haskey(rd.records, key)
# Base.get(rd::RecordDatabase{K, V}, key::K, default::V) where {K, V} = get(rd.records, key, default)
# Base.get(rd::RecordDatabase{K, V}, key::K) where {K, V} = get(rd.records, key)
# function Base.delete!(rd::RecordDatabase{K, V}, key::K) where {K, V}
#     delete!(rd.records, key)
#     rd
# end

function check_uuid(x::RecordDatabase, y::Record)
    id_y = id(y)
    if Base.haskey(x, id_y)
        id_y_string = string(id_y)
        throw(ArgumentError("The ID $id_y_string is already in the data base."))
    end
end

function check_uuid(x::RecordDatabase, y::RecordDatabase)

    keys_x = Set(keys(x))
    keys_y = Set(keys(y))

    duplicated_keys = intersect(keys_x, keys_y)

    if length(duplicated_keys) > 0
        duplicated_keys_str = join(duplicated_keys)
        throw(ArgumentError("Duplicated ID(s): $duplicated_keys_str"))
    end

end
check_uuid(x::RecordDatabase, y::Pair) = check_uuid(x, y.second)


function check_url(x::RecordDatabase, y::Record)
    url_y = url(y)
    if url_y in url(x)
        throw(ArgumentError("The URL $url_y is already in the data base."))
    end
end
function check_url(x::RecordDatabase, y::RecordDatabase)
    duplicated_urls = intersect(url(x), url(y))

    if length(duplicated_urls) > 0
        throw(ArgumentError("Duplicated URL(s): $(join(duplicated_urls))"))
    end
end
check_url(x::RecordDatabase, y::Pair) = check_url(x,y.second)