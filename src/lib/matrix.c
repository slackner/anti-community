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

double *multiply_graph_vector(const struct graph *g, double *vector)
{
    struct link *link;
    double *result;
    double value;
    uint32_t i;

    if (!(result = xmalloc(sizeof(*result) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
    {
        value = 0.0;
        GRAPH_FOR_EACH_LINK(&g->nodes[i], link)
        {
            value += link->weight * vector[link->index];
        }
        result[i] = value;
    }

    return result;
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

struct graph *add_graph_elementwise(const struct graph *g, const struct graph *h)
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
            if (link1 && link2)
                _set_edge(dst, i, link1->index, link1->weight + link2->weight);
            else if (link1)
                _set_edge(dst, i, link1->index, link1->weight);
            else if (link2)
                _set_edge(dst, i, link2->index, link2->weight);
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

double *graph_power_iteration(const struct graph *g, uint32_t num_iterations, double *eigenvalue_out)
{
    double *vector;
    double *temp;
    double norm = 0.0;
    uint32_t i, j;

    if (!num_iterations)
        num_iterations = 100;

    if (!(vector = xmalloc(sizeof(*vector) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
        vector[i] = random_float();

    for (i = 0; i < num_iterations; i++)
    {
        temp = multiply_graph_vector(g, vector);
        free(vector);
        vector = temp;

        norm = 0.0;
        for (j = 0; j < g->num_nodes; j++)
            norm += vector[j] * vector[j];

        norm = sqrt(norm);
        for (j = 0; j < g->num_nodes; j++)
            vector[j] /= norm;
    }

    if (eigenvalue_out)
    {
        temp = multiply_graph_vector(g, vector);
        norm = 0.0;
        for (i = 0; i < g->num_nodes; i++)
            norm += vector[i] * temp[i];
        free(temp);

        *eigenvalue_out = norm;
    }

    return vector;
}
