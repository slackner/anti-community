/*
 * Minimal graph library
 * Adjacency functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include "graph.h"
#include "internal.h"

static int _sort_links_by_index(const void *a, const void *b)
{
    const struct link *la = a, *lb = b;
    return COMPARE(la->index, lb->index);
}

static void adj_init(struct adjacency *adj)
{
    if (adj->num_links)
    {
        qsort(adj->links, adj->num_links, sizeof(adj->links[0]), _sort_links_by_index);
        assert(adj->links[0].index <= adj->links[adj->num_links - 1].index);
    }
}

static struct link *adj_get(struct adjacency *adj, uint32_t end, int allocate)
{
    struct link *link;
    uint32_t insert = 0;

    if (adj->num_links)
    {
        /* Note that we have to use signed numbers here. */
        int32_t min = 0;
        int32_t max = adj->num_links - 1;
        int32_t i;

        if (adj->hint < adj->num_links)
        {
            i = adj->hint;
            link = &adj->links[i];
            if (end < link->index) max = i - 1;
            else if (end > link->index) min = i + 1;
            else
            {
                adj->hint++;
                return link;
            }
        }

        while (min <= max)
        {
            i = (min + max) / 2;
            link = &adj->links[i];
            if (end < link->index) max = i - 1;
            else if (end > link->index) min = i + 1;
            else
            {
                adj->hint = i + 1;
                return link;
            }
        }

        insert = min;
    }

    if (!allocate) return NULL;
    if (!adjacency_reserve(adj, 1)) return NULL;

    link = &adj->links[insert];
    memmove(&link[1], link, (char *)&adj->links[adj->num_links] - (char *)link);
    adj->num_links++;

    link->index = end;
    link->weight = 0.0;
    return link;
}

static void adj_del(struct adjacency *adj, struct link *link)
{
    if (adj->hint < adj->num_links && link <= &adj->links[adj->hint])
    {
        if (link < &adj->links[adj->hint]) adj->hint--;
        else adj->hint = ~0U;
    }
    memmove(link, &link[1], (char *)&adj->links[adj->num_links] - (char *)&link[1]);
    adj->num_links--;
}

const struct adjacency_ops adjacency_ops_sorted =
{
    adj_init,
    adj_get,
    adj_del,
};
