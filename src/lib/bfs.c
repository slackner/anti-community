/*
 * Minimal graph library
 * Breath-first search functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include "graph.h"
#include "internal.h"

static int _sort_bfs_entry_by_weight(const void *a, const void *b, void *userdata)
{
    const struct bfs_entry *ba = a, *bb = b;
    return COMPARE(ba->weight, bb->weight);
}

static int _sort_bfs_entry_by_count(const void *a, const void *b, void *userdata)
{
    const struct bfs_entry *ba = a, *bb = b;
    return COMPARE(ba->count, bb->count);
}

int graph_bfs(const struct graph *g, uint32_t start, int use_weights,
              int (*callback)(const struct graph *g, struct bfs_entry *entry, void *userdata),
              void *userdata)
{
    struct bfs_entry entry, new_entry;
    struct minheap *queue;
    uint32_t i, index;
    struct link *link;
    uint8_t *visited;
    int found = 0;

    if (start >= g->num_nodes)
        return 0;

    if (!(visited = xmalloc(sizeof(*visited) * g->num_nodes)))
        return 0;

    if (!(queue = alloc_minheap(sizeof(struct bfs_entry), use_weights ?
          _sort_bfs_entry_by_weight : _sort_bfs_entry_by_count, NULL)))
    {
        free(visited);
        return 0;
    }

    for (i = 0; i < g->num_nodes; i++)
        visited[i] = 0;

    new_entry.weight = 0.0;
    new_entry.count = 0;
    new_entry.from = ~0U;
    new_entry.to = start;
    if (!minheap_push(queue, &new_entry))
        goto error;

    while (minheap_pop(queue, &entry))
    {
        index = entry.to;
        if (visited[index]) continue;
        if (callback(g, &entry, userdata))
        {
            found = 1;
            break;
        }
        visited[index] = 1;

        GRAPH_FOR_EACH_LINK(&g->nodes[index], link)
        {
            if (visited[link->index]) continue;

            new_entry.weight = entry.weight + link->weight;
            new_entry.count = entry.count + 1;
            new_entry.from = index;
            new_entry.to = link->index;
            if (!minheap_push(queue, &new_entry))
                goto error;
        }
    }

error:
    free_minheap(queue);
    free(visited);
    return found;
}

struct bfs_distance_context
{
    uint32_t end;
    struct bfs_entry entry;
};

static int _bfs_distance_callback(const struct graph *g, struct bfs_entry *entry, void *userdata)
{
    struct bfs_distance_context *context = userdata;
    if (entry->to != context->end) return 0;
    context->entry = *entry;
    return 1;
}

uint32_t graph_get_distance_count(const struct graph *g, uint32_t start, uint32_t end)
{
    struct bfs_distance_context context;
    context.end = end;

    if (graph_bfs(g, start, 0, _bfs_distance_callback, &context))
        return context.entry.count;

    return ~0U;
}

double graph_get_distance_weight(const struct graph *g, uint32_t start, uint32_t end)
{
    struct bfs_distance_context context;
    context.end = end;

    if (graph_bfs(g, start, 1, _bfs_distance_callback, &context))
        return context.entry.weight;

    return INFINITY;
}

struct bfs_all_distances_count_context
{
    uint32_t *counts;
    uint32_t max_count;
};

static int _bfs_all_distances_count_callback(const struct graph *g, struct bfs_entry *entry, void *userdata)
{
    struct bfs_all_distances_count_context *context = userdata;
    if (entry->count > context->max_count) return 1;
    context->counts[entry->to] = entry->count;
    return 0;
}

uint32_t *graph_get_all_distances_count(const struct graph *g, uint32_t start, uint32_t max_count)
{
    struct bfs_all_distances_count_context context;
    uint32_t i;

    context.max_count = max_count;
    if (!(context.counts = xmalloc(sizeof(*context.counts) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++) context.counts[i] = ~0U;
    graph_bfs(g, start, 0, _bfs_all_distances_count_callback, &context);

    return context.counts;
}

struct bfs_all_distances_weight_context
{
    double *weights;
    double max_weight;
};

static int _bfs_all_distances_weight_callback(const struct graph *g, struct bfs_entry *entry, void *userdata)
{
    struct bfs_all_distances_weight_context *context = userdata;
    if (entry->weight > context->max_weight) return 1;
    context->weights[entry->to] = entry->weight;
    return 0;
}

double *graph_get_all_distances_weight(const struct graph *g, uint32_t start, double max_weight)
{
    struct bfs_all_distances_weight_context context;
    uint32_t i;

    context.max_weight = max_weight;
    if (!(context.weights = xmalloc(sizeof(*context.weights) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++) context.weights[i] = INFINITY;
    graph_bfs(g, start, 1, _bfs_all_distances_weight_callback, &context);

    return context.weights;
}

struct bfs_all_distances_graph_context
{
    int use_weights;
    uint32_t start;
    struct graph *distances;
};

static int _bfs_all_distances_graph_callback(const struct graph *g, struct bfs_entry *entry, void *userdata)
{
    struct bfs_all_distances_graph_context *context = userdata;
    if (entry->to == context->start) return 0;  /* skip empty diagonal */
    _add_edge(context->distances, context->start, entry->to, context->use_weights ? entry->weight : entry->count);
    return 0;
}

struct graph *graph_get_all_distances_graph(const struct graph *g, int use_weights)
{
    struct bfs_all_distances_graph_context context;
    uint32_t i;

    context.use_weights = use_weights;
    if (!(context.distances = alloc_graph(g->num_nodes, GRAPH_FLAGS_DIRECTED)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
    {
        context.start = i;
        graph_bfs(g, i, use_weights, _bfs_all_distances_graph_callback, &context);
    }

    return context.distances;
}

struct bfs_components_context
{
    uint32_t *components;
    uint32_t identifier;
};

static int _bfs_components_callback(const struct graph *g, struct bfs_entry *entry, void *userdata)
{
    struct bfs_components_context *context = userdata;
    assert(context->components[entry->to] == ~0U);
    context->components[entry->to] = context->identifier;
    return 0;
}

uint32_t *graph_get_connected_components(const struct graph *g)
{
    struct bfs_components_context context;
    uint32_t i;

    assert(!(g->flags & GRAPH_FLAGS_DIRECTED));

    if (!(context.components = xmalloc(sizeof(*context.components) * g->num_nodes)))
        return NULL;

    context.identifier = 0;
    for (i = 0; i < g->num_nodes; i++) context.components[i] = ~0U;

    for (i = 0; i < g->num_nodes; i++)
    {
        if (context.components[i] != ~0U) continue;
        graph_bfs(g, i, 0, _bfs_components_callback, &context);
        context.identifier++;
    }

    return context.components;
}
