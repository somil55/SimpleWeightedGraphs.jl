__precompile__(true)
module SimpleWeightedGraphs

using LightGraphs

import Base:
    convert, eltype, show, ==, Pair, Tuple, copy, length, start, next, done, issubset, zero

import LightGraphs:
    _NI, _insert_and_dedup!, AbstractGraph, AbstractEdge, AbstractEdgeIter,
    src, dst, edgetype, nv, ne, vertices, edges, is_directed,
    add_vertex!, add_edge!, rem_vertex!, rem_edge!,
    has_vertex, has_edge, in_neighbors, out_neighbors,

    indegree, outdegree, degree, has_self_loops, num_self_loops,

    add_vertices!, adjacency_matrix, weights

import LightGraphs.SimpleGraphs:
    AbstractSimpleEdge, SimpleEdge, fadj, eltype, src, dst, show,
    Pair, Tuple, reverse

export

    SimpleWeightedGraph,
    # SimpleWeightedDiGraph,
    weight,
    weighttype,
    Graph
    # DiGraph,
    # WGraph,
    # WDiGraph

"""
    AbstractSimpleWeightedGraph

An abstract type representing a simple graph structure.
AbstractSimpleWeightedGraphs must have the following elements:
- weightmx::AbstractSparseMatrix{Real}
"""
abstract type AbstractSimpleWeightedGraph <: AbstractGraph end

function show(io::IO, g::AbstractSimpleWeightedGraph)
    if is_directed(g)
        dir = "directed"
    else
        dir = "undirected"
    end
    if nv(g) == 0
        print(io, "empty $dir simple $(eltype(g)) graph with $(weighttype(g)) weights")
    else
        print(io, "{$(nv(g)), $(ne(g))} $dir simple $(eltype(g)) graph with $(weighttype(g)) weights")
    end
end


### INTERFACE
nv(g::AbstractSimpleWeightedGraph) = eltype(g)(length(fadj(g)))
vertices(g::AbstractSimpleWeightedGraph) = one(eltype(g)):nv(g)
weighttype(g::AbstractSimpleWeightedGraph) = eltype(g.weight__.vals)
fadj(g::AbstractSimpleWeightedGraph) = g.fadjlist
fadj(g::AbstractSimpleWeightedGraph, v::Integer) = g.fadjlist[v]


# handles single-argument edge constructors such as pairs and tuples
has_edge(g::AbstractSimpleWeightedGraph, x::Union{Pair,Tuple}) = has_edge(g, edgetype(g)(x))
has_edge(g::AbstractSimpleWeightedGraph, x::Union{Pair,Tuple}, y::Real) = has_edge(g, edgetype(g)(x), y;check_weight=true)
add_edge!(g::AbstractSimpleWeightedGraph, x::Union{Pair,Tuple}) = add_edge!(g, edgetype(g)(x))
add_edge!(g::AbstractSimpleWeightedGraph, x::Union{Pair,Tuple}, y::Real) = add_edge!(g, edgetype(g)(x), y)

# handles two-argument edge constructors like src,dst
has_edge(g::AbstractSimpleWeightedGraph, x::Integer, y::Integer) = has_edge(g, edgetype(g)(x, y))
has_edge(g::AbstractSimpleWeightedGraph, x::Integer, y::Integer, z::Real) = has_edge(g, edgetype(g)(x, y), z;check_weight=true)
add_edge!(g::AbstractSimpleWeightedGraph, x::Integer, y::Integer) = add_edge!(g, edgetype(g)(x, y))
add_edge!(g::AbstractSimpleWeightedGraph, x::Integer, y::Integer, z::Real) = add_edge!(g, edgetype(g)(x, y), z)


in_neighbors(g::AbstractSimpleWeightedGraph, v::Integer) = badj(g,v)
out_neighbors(g::AbstractSimpleWeightedGraph, v::Integer) = fadj(g,v)

has_vertex(g::AbstractSimpleWeightedGraph, v::Integer) = v in vertices(g)

function rem_edge!(g::AbstractSimpleWeightedGraph, u::Integer, v::Integer)
    T = eltype(g)
    rem_edge!(g, edgetype(g)(T(u), T(v)))
end

@doc_str """
    rem_vertex!(g, v)

Remove the vertex `v` from graph `g`. Return false if removal fails
(e.g., if vertex is not in the graph); true otherwise.

### Performance
Time complexity is ``\\mathcal{O}(k^2)``, where ``k`` is the max of the degrees
of vertex ``v`` and vertex ``|V|``.

### Implementation Notes
This operation has to be performed carefully if one keeps external
data structures indexed by edges or vertices in the graph, since
internally the removal is performed swapping the vertices `v`  and ``|V|``,
and removing the last vertex ``|V|`` from the graph. After removal the
vertices in `g` will be indexed by ``1:|V|-1``.
"""
function rem_vertex!(g::AbstractSimpleWeightedGraph, v::Integer)
    v in vertices(g) || return false
    n = nv(g)

    # remove the in_edges from v
    srcs = copy(in_neighbors(g, v))
    for s in srcs
        rem_edge!(g, edgetype(g)(s, v))
    end
    # remove the in_edges from the last vertex
    neigs = copy(in_neighbors(g, n))
    for s in neigs
        rem_edge!(g, edgetype(g)(s, n))
    end
    if v != n
        # add the edges from n back to v
        for s in neigs
            add_edge!(g, edgetype(g)(s, v))
        end
    end

    if is_directed(g)
        # remove the out_edges from v
        dsts = copy(out_neighbors(g, v))
        for d in dsts
            rem_edge!(g, edgetype(g)(v, d))
        end
        # remove the out_edges from the last vertex
        neigs = copy(out_neighbors(g, n))
        for d in neigs
            rem_edge!(g, edgetype(g)(n, d))
        end
        if v != n
            # add the out_edges back to v
            for d in neigs
                add_edge!(g, edgetype(g)(v, d))
            end
        end
    end

    pop!(g.fadjlist)
    if is_directed(g)
        pop!(g.badjlist)
    end
    return true
end


# function issubset(g::T, h::T) where T<:AbstractSimpleWeightedGraph
#     (gmin, gmax) = extrema(vertices(g))
#     (hmin, hmax) = extrema(vertices(h))
#     return (hmin <= gmin <= gmax <= hmax) && issubset(edges(g), edges(h))
# end

# @doc_str """
#     rem_vertex!(g::AbstractSimpleWeightedGraph, v)
#
# Remove the vertex `v` from graph `g`. Return false if removal fails
# (e.g., if vertex is not in the graph); true otherwise.
#
# ### Implementation Notes
# This operation has to be performed carefully if one keeps external
# data structures indexed by edges or vertices in the graph, since
# internally the removal results in all vertices with indices greater than `v`
# being shifted down one.
# """
# function rem_vertex!(g::AbstractSimpleWeightedGraph, v::Integer)
#     warn("Note: removing vertices from this graph type is not performant.", once=true)
#     v in vertices(g) || return false
#     n = nv(g)
#
#     newweights = g.weights[1:nv(g) .!= v, :]
#     newweights = newweights[:, 1:nv(g) .!= v]
#
#     g.weights = newweights
#     return true
# end

# function out_neighbors(g::AbstractSimpleWeightedGraph)
#     mat = g.weights
#     return [mat.rowval[mat.colptr[i]:mat.colptr[i+1]-1] for i in 1:nv(g)]
# end
#
# function out_neighbors(g::AbstractSimpleWeightedGraph, v::Integer)
#     mat = g.weights
#     return mat.rowval[mat.colptr[v]:mat.colptr[v+1]-1]
# end
#
# zero(g::T) where T<:AbstractSimpleWeightedGraph = T()
#
# # TODO: manipulte SparseMatrixCSC directly
# add_vertex!(g::AbstractSimpleWeightedGraph) = add_vertices!(g, 1)
#
# copy(g::T) where T <: AbstractSimpleWeightedGraph =  T(copy(g.weights))


# const SimpleWeightedGraphEdge = SimpleWeightedEdge
# const SimpleWeightedDiGraphEdge = SimpleWeightedEdge
# include("simpleweighteddigraph.jl")
include("simpleweightedgraph.jl")
# include("overrides.jl")


# const Graph = SimpleWeightedGraph
# const WGraph = SimpleWeightedGraph
# const WDiGraph = SimpleWeightedDiGraph

end # module
