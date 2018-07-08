/*
 * Minimal graph library
 * Metric functions.
 *
 * Copyright (c) 2018 Sebastian Lackner
 */

#include "graph.h"
#include "internal.h"

double *graph_closeness_centrality(const struct graph *g)
{
    double *centrality;
    uint32_t *counts;
    double distances;
    uint32_t i, j;

    if (!(centrality = xmalloc(sizeof(*centrality) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
    {
        if (!(counts = graph_get_all_distances_count(g, i, ~0U)))
        {
            free(centrality);
            return NULL;
        }

        distances = 0.0;
        for (j = 0; j < g->num_nodes; j++)
        {
            if (i == j) continue;
            distances += counts[j];
        }

        centrality[i] = (double)(g->num_nodes - 1) / distances;
        free(counts);
    }

    return centrality;
}

double *graph_harmonic_centrality(const struct graph *g, uint32_t max_count)
{
    double *centrality;
    uint32_t *counts;
    double distances;
    uint32_t i, j;

    if (!(centrality = xmalloc(sizeof(*centrality) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
    {
        if (!(counts = graph_get_all_distances_count(g, i, max_count)))
        {
            free(centrality);
            return NULL;
        }

        distances = 0.0;
        for (j = 0; j < g->num_nodes; j++)
        {
            if (i == j) continue;
            distances += 1.0 / counts[j];
        }

        centrality[i] = distances / (double)(g->num_nodes - 1);
        free(counts);
    }

    return centrality;
}

double *graph_harmonic_centrality_fast(const struct graph *g, uint32_t max_count)
{
    struct link *link1, *link2;
    double *centrality;
    double distances;
    uint8_t *visited;
    uint32_t i, j;

    assert(max_count <= 3);

    if (!(visited = xmalloc(sizeof(*visited) * g->num_nodes)))
        return NULL;

    if (!(centrality = xmalloc(sizeof(*centrality) * g->num_nodes)))
    {
        free(visited);
        return NULL;
    }

    for (i = 0; i < g->num_nodes; i++)
    {
        for (j = 0; j < g->num_nodes; j++)
            visited[j] = 0;

        visited[i] = 1;
        distances = 0;

        if (max_count < 1) goto done;
        GRAPH_FOR_EACH_LINK(&g->nodes[i], link1)
        {
            if (link1->index == i) continue;
            visited[link1->index] = 1;
            distances += 1.0;
        }

        if (max_count < 2) goto done;
        GRAPH_FOR_EACH_LINK(&g->nodes[i], link1)
        {
            if (link1->index == i) continue;
            GRAPH_FOR_EACH_LINK(&g->nodes[link1->index], link2)
            {
                if (visited[link2->index]) continue;
                visited[link2->index] = 2;
                distances += (1.0 / 2.0);
            }
        }

        if (max_count < 3) goto done;
        for (j = 0; j < g->num_nodes; j++)
        {
            if (visited[j] != 2) continue;
            GRAPH_FOR_EACH_LINK(&g->nodes[j], link2)
            {
                if (visited[link2->index]) continue;
                visited[link2->index] = 3;
                distances += (1.0 / 3.0);
            }
        }

done:
        centrality[i] = distances / (double)(g->num_nodes - 1);
    }

    free(visited);
    return centrality;
}

double graph_transitivity(const struct graph *g)
{
    struct link *link1, *link2, *link3;
    uint64_t triangles = 0;
    uint64_t triples = 0;
    uint32_t i, count;

    assert(!(g->flags & GRAPH_FLAGS_DIRECTED));

    for (i = 0; i < g->num_nodes; i++)
    {
        count = 0;
        GRAPH_FOR_EACH_LINK(&g->nodes[i], link1)
        {
            if (link1->index == i) continue;
            GRAPH_FOR_EACH_LINK2(&g->nodes[link1->index], link2, &g->nodes[i], link3)
            {
                if (!link2 || !link3) continue;
                if (link2->index == i) continue;
                if (link2->index == link1->index) continue;
                /* triangle between (i, link1->index, link2->index) */
                triangles++;
            }
            count++;
        }
        triples += (uint64_t)count * (count - 1);
    }

    if (!triangles)
        return 0.0;

    return (double)triangles / triples;
}

double *graph_clustering_coefficient(const struct graph *g, int use_weights)
{
    struct link *link1, *link2, *link3;
    uint64_t triangles, triples;
    uint32_t i, count;
    double *clustering;

    assert(!(g->flags & GRAPH_FLAGS_DIRECTED));
    assert(!use_weights);  /* FIXME: not implemented yet */

    if (!(clustering = xmalloc(sizeof(*clustering) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
    {
        triangles = 0;
        count = 0;
        GRAPH_FOR_EACH_LINK(&g->nodes[i], link1)
        {
            if (link1->index == i) continue;
            GRAPH_FOR_EACH_LINK2(&g->nodes[link1->index], link2, &g->nodes[i], link3)
            {
                if (!link2 || !link3) continue;
                if (link2->index == i) continue;
                if (link2->index == link1->index) continue;
                /* triangle between (i, link1->index, link2->index) */
                triangles++;
            }
            count++;
        }
        triples = (uint64_t)count * (count - 1);
        clustering[i] = (triangles != 0) ? ((double)triangles / triples) : 0.0;
    }

    return clustering;
}

double graph_reciprocity(const struct graph *g)
{
    uint64_t count_total = 0;
    uint64_t count_bidir = 0;
    struct link *link;
    uint32_t i;

    assert(g->flags & GRAPH_FLAGS_DIRECTED);

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        if (_get_link(&g->nodes[link->index], i, 0)) count_bidir++;
        count_total++;
    }

    return (double)count_bidir / count_total;
}
