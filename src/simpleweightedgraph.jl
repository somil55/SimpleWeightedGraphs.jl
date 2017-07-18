"""
    SimpleWeightedGraph{T, U}

A type representing an undirected graph with weights of type `U`.
"""
mutable struct SimpleWeightedGraph{T<:Integer, U<:Real} <: AbstractSimpleWeightedGraph
  ne::Int
  fadjlist::Vector{Vector{T}} # [src]: (dst, dst, dst)
  weight__::Dict{LightGraphs.SimpleGraphs.SimpleEdge{T},U}
end

ne(g::SimpleWeightedGraph) = g.ne


# Graph{UInt8}(6), Graph{Int16}(7), Graph{UInt8}()
function (::Type{SimpleWeightedGraph{T, U}})(n::Integer = 0) where T<:Integer where U<:Real
    fadjlist = Vector{Vector{T}}()
    sizehint!(fadjlist, n)
    for _ = one(T):n
        push!(fadjlist, Vector{T}())
    end
    weight__ = Dict{LightGraphs.SimpleGraphs.SimpleEdge{T},U}()
    return SimpleWeightedGraph{T, U}(0, fadjlist, weight__)
end


# Graph()
SimpleWeightedGraph() = SimpleWeightedGraph{Int, Float64}()

# Graph(6), Graph(0x5)
SimpleWeightedGraph(n::T) where T<:Integer = SimpleWeightedGraph{T, Float64}(n)

# Graph(UInt8)
SimpleWeightedGraph(::Type{T}) where T<:Integer = SimpleWeightedGraph{T, Float64}(zero(T))

# Graph(UInt8, Float32)
SimpleWeightedGraph(::Type{T}, ::Type{U}) where T<:Integer where U<:Real = SimpleWeightedGraph{U, T}(zero(T))

eltype(g::SimpleWeightedGraph{T,U}) where T<:Integer where U<:Real = T
edgetype(::SimpleWeightedGraph{T,U}) where T<:Integer where U<:Real = SimpleEdge{T}

# edges(g::SimpleWeightedGraph) = (SimpleWeightedEdge(x[1], x[2], x[3]) for x in zip(findnz(triu(g.weights))...))
weights(g::SimpleWeightedGraph) = g.weight__

function has_edge(
  g::SimpleWeightedGraph,
  e::SimpleEdge,
  w::Real=0;
  check_weight=false)

    u, v = Tuple(e)
    u > nv(g) || v > nv(g) && return false
    if degree(g,u) > degree(g,v)
        u, v = v, u
    end
    edge_found = length(searchsorted(fadj(g,u), v)) > 0
    if edge_found && check_weight
      edge_found = (g.weight__[e] == w)
    end

    return edge_found
end

function add_edge!(
  g::SimpleWeightedGraph,
  e::LightGraphs.SimpleGraphs.SimpleEdge,
  w::Real=1.0)
    T = eltype(g)
    U = weighttype(g)
    s, d = T.(Tuple(e))
    (s in vertices(g) && d in vertices(g)) || return false
    inserted = _insert_and_dedup!(g.fadjlist[s], d)
    if inserted
        g.ne += 1
    end
    if s != d
        inserted = _insert_and_dedup!(g.fadjlist[d], s)
    end
    g.weight__[e] = U.(w)
    g.weight__[reverse(e)] = U.(w)
    return inserted
end

function rem_edge!(g::AbstractSimpleWeightedGraph, e::SimpleEdge)
  U = weighttype(g)
  i = searchsorted(g.fadjlist[src(e)], dst(e))
  length(i) > 0 || return false   # edge not in graph
  i = i[1]
  deleteat!(g.fadjlist[src(e)], i)
  if src(e) != dst(e)     # not a self loop
    i = searchsorted(g.fadjlist[dst(e)], src(e))[1]
    deleteat!(g.fadjlist[dst(e)], i)
  end
  delete!(g.weight__, e)
  delete!(g.weight__, reverse(e))
  g.ne -= 1

  return true
end

"""
    is_directed(g)

Return `true` if `g` is a directed graph.
"""
is_directed(::Type{SimpleWeightedGraph}) = false
is_directed(::Type{SimpleWeightedGraph{T, U}}) where T where U = false
is_directed(g::SimpleWeightedGraph) = false


"""
    add_vertex!(g)

Add a new vertex to the graph `g`. Return `true` if addition was successful.
"""
function add_vertex!(g::SimpleWeightedGraph)
    T = eltype(g)
    (nv(g) + one(T) <= nv(g)) && return false       # test for overflow
    push!(g.fadjlist, Vector{T}())
    return true
end

"""
    badj(g::SimpleGraph[, v::Integer])

Return the backwards adjacency list of a graph. If `v` is specified,
return only the adjacency list for that vertex.

###Implementation Notes
Returns a reference, not a copy. Do not modify result.
"""
badj(g::SimpleWeightedGraph) = fadj(g)
badj(g::SimpleWeightedGraph, v::Integer) = fadj(g, v)

# # Graph(SimpleGraph)
# SimpleWeightedGraph(g::LightGraphs.SimpleGraphs.SimpleGraph{T}, ::Type{U}=Float64) where T <: Integer where U <: Real =
#     SimpleWeightedGraph{T, U}(adjacency_matrix(g))
#
# # Graph(SimpleDiGraph)
# SimpleWeightedGraph(g::LightGraphs.SimpleGraphs.SimpleDiGraph{T}, ::Type{U}=Float64) where T <: Integer where U <: Real =
#     SimpleWeightedGraph{T, U}(adjacency_matrix(LightGraphs.SimpleGraphs.SimpleGraph(g)))
#
# # Graph(SimpleGraph, defaultweight)
# SimpleWeightedGraph(g::LightGraphs.SimpleGraphs.SimpleGraph{T}, x::U) where T <: Integer where U <: Real =
#     SimpleWeightedGraph{T, U}(x.*adjacency_matrix(g, U))
#
# # Graph(SimpleDiGraph, defaultweight)
# SimpleWeightedGraph(g::LightGraphs.SimpleGraphs.SimpleDiGraph{T}, x::U) where T <: Integer where U <: Real =
#     SimpleWeightedGraph{T, U}(x.*adjacency_matrix(LightGraphs.SimpleGraphs.SimpleGraph(g), U))
#
# # DiGraph(srcs, dsts, weights)
# SimpleWeightedGraph(i::AbstractVector{T}, j::AbstractVector{T}, v::AbstractVector{U}) where T<:Integer where U<:Real =
#     SimpleWeightedGraph{T, U}(sparse(i,j,v))
# Graph{UInt8}(adjmx)
# function (::Type{SimpleWeightedGraph{T, U}})(adjmx::AbstractMatrix) where T<:Integer where U <: Real
#     dima,dimb = size(adjmx)
#     isequal(dima,dimb) || error("Adjacency / distance matrices must be square")
#     issymmetric(adjmx) || error("Adjacency / distance matrices must be symmetric")
#     g = SimpleWeightedGraph(U.(spones(adjmx)))
# end

# converts Graph{Int} to Graph{Int32}
# function (::Type{SimpleWeightedGraph{T, U}})(g::SimpleWeightedGraph) where T<:Integer where U<:Real
#     h_fadj = [Vector{T}(x) for x in fadj(g)]
#     return SimpleGraph(ne(g), h_fadj)
# end


# Graph(adjmx)
# function SimpleWeightedGraph(adjmx::AbstractMatrix)
#     dima,dimb = size(adjmx)
#     isequal(dima,dimb) || error("Adjacency / distance matrices must be square")
#     issymmetric(adjmx) || error("Adjacency / distance matrices must be symmetric")
#     SimpleWeightedGraph{Int, eltype(adjmx)}(adjmx')
# end

# Graph(digraph). Weights will be added.

# SimpleWeightedGraph(g::SimpleWeightedDiGraph) = SimpleWeightedGraph(g.weights .+ g.weights')

# edgetype(::SimpleWeightedGraph{T, U}) where T<:Integer where U<:Real= SimpleWeightedGraphEdge{T,U}



#
#
# ==(g::SimpleWeightedGraph, h::SimpleWeightedGraph) = g.weights == h.weights
#
#
