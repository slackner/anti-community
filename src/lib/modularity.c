/*
 * Minimal graph library
 * Modularity / Anti-Modularity functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#define _GNU_SOURCE
#include "graph.h"
#include "internal.h"

/*
 * Implementation based on:
 *
 * - L. Chen, Q. Yu, and B. Chen. 2014. Anti-modularity and anti-community
 *   detecting in complex networks. Inf. Sci. 275 (2014), 293–313.
 *   https://doi.org/10.1016/j.ins.2014.02.040
 *
 * - M. E. J. Newman and M. Girvan. 2004. Finding and evaluating community
 *   structure in networks. Phys. Rev. E 69, 2 (2004), 15.
 *   https://doi.org/10.1103/PhysRevE.69.026113
 *
 */

static double degree_product(double *norm_out, const struct graph *g, const uint32_t *labels)
{
    double *weights_out, *weights_in;
    uint32_t *lookup;
    double product = 0.0;
    uint32_t i;

    if (!(lookup = xmalloc(sizeof(*lookup) * g->num_nodes)))
        return NAN;

    if (!(weights_out = graph_out_weights(g)))
    {
        free(lookup);
        return NAN;
    }

    if (!(g->flags & GRAPH_FLAGS_DIRECTED)) weights_in = weights_out;
    else if (!(weights_in = graph_in_weights(g)))
    {
        free(weights_out);
        free(lookup);
        return NAN;
    }

    for (i = 0; i < g->num_nodes; i++) lookup[i] = i;
    qsort_r(lookup, g->num_nodes, sizeof(lookup[0]), _sort_index_by_label, (void *)labels);

    for (i = 0; i < g->num_nodes;)
    {
        uint32_t label = labels[lookup[i]];
        double out = 0.0, in = 0.0;

        for (; i < g->num_nodes; i++)
        {
            if (labels[lookup[i]] != label) break;
            out += weights_out[lookup[i]];
            in += weights_in[lookup[i]];
        }

        product += out * in;
    }

    if (norm_out)
    {
        double norm = 0.0;
        for (i = 0; i < g->num_nodes; i++)
            norm += weights_out[i];
        *norm_out = norm;
    }

    if (weights_in != weights_out) free(weights_in);
    free(weights_out);
    free(lookup);
    return product;
}

double modularity(const struct graph *g, const uint32_t *labels)
{
    double mod[2] = { 0.0, 0.0 };
    struct link *link;
    uint32_t i;
    double norm = NAN;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        if (labels[i] != labels[link->index]) continue;
        mod[0] += link->weight;
    }

    mod[1] = degree_product(&norm, g, labels);
    return (mod[0] - mod[1] / norm) / norm;
}

double modularity_decode(const struct graph *g, const uint32_t *indices)
{
    uint32_t *labels;
    double value = NAN;

    if ((labels = decode_labels(indices, g->num_nodes)))
    {
        value = modularity(g, labels);
        free(labels);
    }

    return value;
}

double antimodularity(const struct graph *g, const uint32_t *labels)
{
    double norm = g->num_nodes, mod[2] = { 0.0, 0.0 };
    struct link *link1, *link2;
    uint32_t i;

    if (g->flags & GRAPH_FLAGS_DIRECTED)
    {
        fprintf(stderr, "FIXME: Anti-modularity not defined for directed graphs\n");
        return NAN;
    }

    GRAPH_FOR_EACH_EDGE(g, i, link1)
    {
        GRAPH_FOR_EACH_LINK(&g->nodes[link1->index], link2)
        {
            if (labels[i] != labels[link2->index]) continue;
            mod[0] += link1->weight * link2->weight;
        }
    }

    mod[1] = degree_product(NULL, g, labels);
    return (mod[0] - mod[1] / norm) / norm;
}

double antimodularity_decode(const struct graph *g, const uint32_t *indices)
{
    uint32_t *labels;
    double value = NAN;

    if ((labels = decode_labels(indices, g->num_nodes)))
    {
        value = antimodularity(g, labels);
        free(labels);
    }

    return value;
}
