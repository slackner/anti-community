/*
 * Minimal graph library
 * Node merge functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include "graph.h"
#include "internal.h"

void graph_merge_nodes_add(struct graph *g, uint32_t index0, uint32_t index1)
{
    struct adjacency *adj;
    struct link *link;
    uint32_t i;

    if (index0 >= g->num_nodes) return;
    if (index1 >= g->num_nodes) return;
    if (index0 == index1) return;

    adj = &g->nodes[index1];

    if (g->flags & GRAPH_FLAGS_DIRECTED)
    {
        if (adj->links)
        {
            GRAPH_FOR_EACH_LINK(adj, link)
            {
                _add_edge(g, index0, link->index, link->weight);
            }
            free_adjacency(adj);
        }

        for (i = 0; i < g->num_nodes; i++)
        {
            if ((link = _get_link(&g->nodes[i], index1, 0)))
            {
                float weight = link->weight;
                _del_link(&g->nodes[i], link);
                _add_edge(g, i, index0, weight);
            }
        }
    }
    else if (adj->links)
    {
        GRAPH_FOR_EACH_LINK(adj, link)
        {
            if (link->index == index0)
            {
                _add_edge(g, index0, index0, 2.0 * link->weight);
                _del_edge(g, link->index, index1);
            }
            else if (link->index == index1)
            {
                _add_edge(g, index0, index0, link->weight);
            }
            else
            {
                _add_edge(g, index0, link->index, link->weight);
                _add_edge(g, link->index, index0, link->weight);
                _del_edge(g, link->index, index1);
            }
        }
        free_adjacency(adj);
    }
}

void graph_merge_nodes_min(struct graph *g, uint32_t index0, uint32_t index1)
{
    struct adjacency *adj;
    struct link *link;
    uint32_t i;

    if (index0 >= g->num_nodes) return;
    if (index1 >= g->num_nodes) return;
    if (index0 == index1) return;

    adj = &g->nodes[index1];

    if (g->flags & GRAPH_FLAGS_DIRECTED)
    {
        if (adj->links)
        {
            GRAPH_FOR_EACH_LINK(adj, link)
            {
                _min_edge(g, index0, link->index, link->weight);
            }
            free_adjacency(adj);
        }

        for (i = 0; i < g->num_nodes; i++)
        {
            if ((link = _get_link(&g->nodes[i], index1, 0)))
            {
                float weight = link->weight;
                _del_link(&g->nodes[i], link);
                _min_edge(g, i, index0, weight);
            }
        }
    }
    else if (adj->links)
    {
        GRAPH_FOR_EACH_LINK(adj, link)
        {
            if (link->index == index0)
            {
                _min_edge(g, index0, index0, link->weight);
                _del_edge(g, link->index, index1);
            }
            else if (link->index == index1)
            {
                _min_edge(g, index0, index0, link->weight);
            }
            else
            {
                _min_edge(g, index0, link->index, link->weight);
                _min_edge(g, link->index, index0, link->weight);
                _del_edge(g, link->index, index1);
            }
        }
        free_adjacency(adj);
    }
}

void graph_merge_nodes_max(struct graph *g, uint32_t index0, uint32_t index1)
{
    struct adjacency *adj;
    struct link *link;
    uint32_t i;

    if (index0 >= g->num_nodes) return;
    if (index1 >= g->num_nodes) return;
    if (index0 == index1) return;

    adj = &g->nodes[index1];

    if (g->flags & GRAPH_FLAGS_DIRECTED)
    {
        if (adj->links)
        {
            GRAPH_FOR_EACH_LINK(adj, link)
            {
                _max_edge(g, index0, link->index, link->weight);
            }
            free_adjacency(adj);
        }

        for (i = 0; i < g->num_nodes; i++)
        {
            if ((link = _get_link(&g->nodes[i], index1, 0)))
            {
                float weight = link->weight;
                _del_link(&g->nodes[i], link);
                _max_edge(g, i, index0, weight);
            }
        }
    }
    else if (adj->links)
    {
        GRAPH_FOR_EACH_LINK(adj, link)
        {
            if (link->index == index0)
            {
                _max_edge(g, index0, index0, link->weight);
                _del_edge(g, link->index, index1);
            }
            else if (link->index == index1)
            {
                _max_edge(g, index0, index0, link->weight);
            }
            else
            {
                _max_edge(g, index0, link->index, link->weight);
                _max_edge(g, link->index, index0, link->weight);
                _del_edge(g, link->index, index1);
            }
        }
        free_adjacency(adj);
    }
}

void graph_merge_nodes_custom(struct graph *g, uint32_t index0, uint32_t index1,
                              double alpha1, double alpha2, double beta, double gamma)
{
    struct adjacency *adj;
    struct link *link0, *link1;
    float weight, constant;

    assert(!(g->flags & GRAPH_FLAGS_DIRECTED));
    if (index0 >= g->num_nodes) return;
    if (index1 >= g->num_nodes) return;
    if (index0 == index1) return;

    adj = &g->nodes[index1];

    if (adj->links)
    {
        constant = beta * graph_get_edge(g, index0, index1);

        GRAPH_FOR_EACH_LINK(adj, link1)
        {
            if (link1->index == index1) continue;
            if (link1->index != index0 && (link0 = _get_link(&g->nodes[index0], link1->index, 1)))
            {
                weight = alpha1 * link0->weight + alpha2 * link1->weight +
                         constant + gamma * abs(link0->weight - link1->weight);

                link0->weight = weight;
                _set_edge(g, link1->index, index0, weight);
            }
            _del_edge(g, link1->index, index1);
        }

        _del_edge(g, index0, index0);
        free_adjacency(adj);
    }
}
