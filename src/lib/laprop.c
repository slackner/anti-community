/*
 * Minimal graph library
 * Label propagation algorithm for anti-community detection.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include "graph.h"

static int label_compat(struct adjacency *adj, uint32_t *labels, uint32_t new_label)
{
    struct link *link;

    GRAPH_FOR_EACH_LINK(adj, link)
    {
        if (labels[link->index] == new_label) return 0;
    }

    return 1;
}

uint32_t *laprop(struct graph *g)
{
    uint32_t *labels = NULL;
    uint32_t i, j;
    uint32_t num_iter;
    int changed = 1;

    assert(!(g->flags & GRAPH_FLAGS_DIRECTED));

    labels = alloc_labels(g->num_nodes);

    for (num_iter = 0; changed; num_iter++)                                     /* O(n) */
    {
        progress("Iteration %u", num_iter);

        changed = 0;
        for (i = 0; i < g->num_nodes; i++)
        {
            GRAPH_FOR_EACH_NONEXISTENT_LINK(g, &g->nodes[i], j)                  /* O(n) */
            {
                if (labels[j] < labels[i] && label_compat(&g->nodes[i], labels, labels[j]))
                {
                    labels[i] = labels[j];
                    changed = 1;
                    break;
                }
            }
        }
    }

    return labels;
}
