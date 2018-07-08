/*
 * Minimal graph library
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#ifndef _INTERNAL_H_
#define _INTERNAL_H_

#include "graph.h"

#define MIN(a, b) \
    ({ __typeof__ (a) __a = (a); \
       __typeof__ (b) __b = (b); \
       __a < __b ? __a : __b; })

#define MAX(a, b) \
    ({ __typeof__ (a) __a = (a); \
       __typeof__ (b) __b = (b); \
       __a > __b ? __a : __b; })

#define COMPARE(a, b) \
    ({ __typeof__ (a) __a = (a); \
       __typeof__ (b) __b = (b); \
       (__a > __b) - (__a < __b); })

#define SWAP(a, b, size)                        \
    do                                          \
    {                                           \
        register size_t __size = (size);        \
        register char *__a = (a), *__b = (b);   \
        while (__size--)                        \
        {                                       \
            char __tmp = *__a;                  \
            *__a++ = *__b;                      \
            *__b++ = __tmp;                     \
        }                                       \
    }                                           \
    while (0)

#define LIKELY(x)   __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)
#define DECL_INTERNAL __attribute__((__visibility__("hidden")))

static inline struct link *_get_link(struct adjacency *adj, uint32_t end, int allocate)
{
    return adj->ops->get(adj, end, allocate);
}

static inline void _del_link(struct adjacency *adj, struct link *link)
{
    adj->ops->del(adj, link);
}

static inline void _set_edge(struct graph *g, uint32_t start, uint32_t end, float weight)
{
    struct link *link;
    assert(start < g->num_nodes && end < g->num_nodes);
    if ((link = _get_link(&g->nodes[start], end, 1))) link->weight = weight;
}

static inline void _add_edge(struct graph *g, uint32_t start, uint32_t end, float weight)
{
    struct link *link;
    assert(start < g->num_nodes && end < g->num_nodes);
    if ((link = _get_link(&g->nodes[start], end, 1))) link->weight += weight;
}

static inline void _min_edge(struct graph *g, uint32_t start, uint32_t end, float weight)
{
    struct link *link;
    assert(start < g->num_nodes && end < g->num_nodes);
    if ((link = _get_link(&g->nodes[start], end, 1))) link->weight = MIN(link->weight, weight);
}

static inline void _max_edge(struct graph *g, uint32_t start, uint32_t end, float weight)
{
    struct link *link;
    assert(start < g->num_nodes && end < g->num_nodes);
    if ((link = _get_link(&g->nodes[start], end, 1))) link->weight = MAX(link->weight, weight);
}

static inline void _add_edge_swap(struct graph *g, uint32_t start, uint32_t end, float weight)
{
    struct link *link;
    assert(start < g->num_nodes && end < g->num_nodes);
    if ((link = _get_link(&g->nodes[end], start, 1))) link->weight += weight;
}

static inline void _del_edge(struct graph *g, uint32_t start, uint32_t end)
{
    struct link *link;
    assert(start < g->num_nodes && end < g->num_nodes);
    if ((link = _get_link(&g->nodes[start], end, 0))) _del_link(&g->nodes[start], link);
}

static inline int _sort_index_by_label(const void *a, const void *b, void *userdata)
{
    const uint32_t *labels = userdata;
    const uint32_t *ia = a, *ib = b;
    return COMPARE(labels[*ia], labels[*ib]);
}

struct sort_index_by_label2_context
{
    const uint32_t *labels1;
    const uint32_t *labels2;
};

static inline int _sort_index_by_label2(const void *a, const void *b, void *userdata)
{
    const struct sort_index_by_label2_context *context = userdata;
    const uint32_t *ia = a, *ib = b;
    int res;

    if ((res = COMPARE(context->labels1[*ia], context->labels1[*ib]))) return res;
    return COMPARE(context->labels2[*ia], context->labels2[*ib]);
}

static inline void init_adjacency(struct adjacency *adj, const struct adjacency_ops *ops)
{
    adj->ops       = ops;
    adj->links     = NULL;
    adj->num_links = 0;
    adj->max_links = 0;
    adj->hint      = ~0U;
}

static inline void free_adjacency(struct adjacency *adj)
{
    free(adj->links);
    adj->links     = NULL;
    adj->num_links = 0;
    adj->max_links = 0;
    adj->hint      = ~0U;
}

static inline int adjacency_reserve(struct adjacency *adj, uint32_t new_links)
{
    struct link *link;
    uint32_t max_links;

    if (!adj->links)
    {
        max_links = MAX(new_links, 2);
        if (!(link = xmalloc(sizeof(*link) * max_links)))
            return 0;

        adj->num_links = 0;
        adj->max_links = max_links;
        adj->links     = link;
    }
    else if (adj->num_links + new_links > adj->max_links)
    {
        max_links = MAX(adj->num_links + new_links, adj->max_links * 2);
        if (!(link = xrealloc(adj->links, sizeof(*link) * max_links)))
            return 0;

        adj->max_links = max_links;
        adj->links     = link;
    }

    return 1;
}

#endif /* _INTERNAL_H_ */
