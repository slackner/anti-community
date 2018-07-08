/*
 * Minimal graph library
 * Matrix functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include "graph.h"
#include "internal.h"

struct graph *multiply_graph_const(const struct graph *g, double constant)
{
    struct graph *dst;
    struct link *link;
    uint32_t i;

    if (!(dst = duplicate_graph(g)))
        return NULL;

    GRAPH_FOR_EACH_EDGE(dst, i, link)
    {
        link->weight *= constant;
    }

    compress_graph_inline(dst);
    return dst;
}

struct graph *multiply_graph_elementwise(const struct graph *g, const struct graph *h)
{
    struct link *link1, *link2;
    struct graph *dst;
    uint32_t i;

    assert(g->num_nodes == h->num_nodes);

    if (!(dst = alloc_graph(g->num_nodes, (g->flags | h->flags) & ~GRAPH_FLAGS_UNSORTED)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
    {
        GRAPH_FOR_EACH_LINK2(&g->nodes[i], link1, &h->nodes[i], link2)
        {
            if (!link1 || !link2) continue;
            _set_edge(dst, i, link1->index, link1->weight * link2->weight);
        }
    }

    compress_graph_inline(dst);
    return dst;
}

double scalar_product_graph(const struct graph *g, const struct graph *h)
{
    struct link *link1, *link2;
    double result = 0.0;
    uint32_t i;

    assert(g->num_nodes == h->num_nodes);

    for (i = 0; i < g->num_nodes; i++)
    {
        GRAPH_FOR_EACH_LINK2(&g->nodes[i], link1, &h->nodes[i], link2)
        {
            if (!link1 || !link2) continue;
            result += link1->weight * link2->weight;
        }
    }

    return result;
}

struct graph *multiply_graph(const struct graph *g, const struct graph *h)
{
    void (*add_edge)(struct graph *, uint32_t, uint32_t, float) = _add_edge;
    struct link *link1, *link2;
    struct graph *dst;
    uint32_t i;

    assert(g->num_nodes == h->num_nodes);

    if (!(dst = alloc_graph(g->num_nodes, (g->flags | h->flags) & ~GRAPH_FLAGS_UNSORTED)))
        return NULL;

    /* Both input graphs are undirected. This means the output graph will be
     * undirected aswell. We can optimize performance by swapping indices. */
    if (!(dst->flags & GRAPH_FLAGS_DIRECTED))
        add_edge = _add_edge_swap;

    GRAPH_FOR_EACH_EDGE(g, i, link1)
    {
        GRAPH_FOR_EACH_LINK(&h->nodes[link1->index], link2)
        {
            add_edge(dst, i, link2->index, link1->weight * link2->weight);
        }
    }

    compress_graph_inline(dst);
    return dst;
}

struct graph *square_graph(const struct graph *g)
{
    return multiply_graph(g, g);
}
