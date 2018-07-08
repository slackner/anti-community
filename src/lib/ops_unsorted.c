/*
 * Minimal graph library
 * Adjacency functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include "graph.h"
#include "internal.h"

static void adj_init(struct adjacency *adj)
{
    /* nothing to do */
}

static struct link *adj_get(struct adjacency *adj, uint32_t end, int allocate)
{
    struct link *link;

    GRAPH_FOR_EACH_LINK(adj, link)
    {
        if (link->index == end) return link;
    }

    if (!allocate) return NULL;
    if (!adjacency_reserve(adj, 1)) return NULL;

    link = &adj->links[adj->num_links++];
    link->index = end;
    link->weight = 0.0;
    return link;
}

static void adj_del(struct adjacency *adj, struct link *link)
{
    *link = adj->links[--adj->num_links];
}

const struct adjacency_ops adjacency_ops_unsorted =
{
    adj_init,
    adj_get,
    adj_del,
};
