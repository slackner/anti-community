/*
 * Minimal graph library
 * Generic functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#define _GNU_SOURCE
#include "graph.h"
#include "internal.h"

void *xmalloc(size_t size)
{
    void *p = malloc(size);
    if (UNLIKELY(!p))
    {
        fprintf(stderr, "out of memory!\n");
        exit(1);
    }
    return p;
}

void *xrealloc(void *ptr, size_t new_size)
{
    void *p = realloc(ptr, new_size);
    if (UNLIKELY(!p))
    {
        fprintf(stderr, "out of memory!\n");
        exit(1);
    }
    return p;
}

struct graph *alloc_graph(uint32_t num_nodes, uint32_t flags)
{
    const struct adjacency_ops *ops = &adjacency_ops_sorted;
    struct graph *g;
    uint32_t i;

    if (flags & GRAPH_FLAGS_UNSORTED)
        ops = &adjacency_ops_unsorted;

    if (!(g = xmalloc(offsetof(struct graph, buffer[num_nodes]))))
        return NULL;

    g->flags        = flags;
    g->num_nodes    = num_nodes;
    g->userdata     = NULL;
    g->nodes        = num_nodes ? g->buffer : NULL;

    for (i = 0; i < num_nodes; i++)
        init_adjacency(&g->nodes[i], ops);

    return g;
}

struct graph *load_graph(const char *filename, uint32_t flags)
{
    struct graph *g = NULL;
    ssize_t read;
    size_t len = 0;
    char *line = NULL;
    FILE *fp;
    int res;

    if (!(fp = fopen(filename, "r")))
        return NULL;

    while ((read = getline(&line, &len, fp)) > 0)
    {
        if (line[read - 1] == '\n') line[read - 1] = 0;
        if (!line[0] || line[0] == '#' || line[0] == ';') continue;
        if (!g)
        {
            uint32_t num_nodes;
            if (sscanf(line, "%u", &num_nodes) < 1) goto error;
            if (!(g = alloc_graph(num_nodes, flags))) goto error;
        }
        else
        {
            uint32_t i, j;
            float weight;
            if ((res = sscanf(line, "%u %u %f", &i, &j, &weight)) < 2) goto error;
            if (i < 1 || i > g->num_nodes) goto error;
            if (j < 1 || j > g->num_nodes) goto error;
            if (res < 3) weight = 1.0;
            graph_set_edge(g, i - 1, j - 1, weight);
        }
    }
    if (!g) goto error;

    fclose(fp);
    free(line);

    compress_graph_inline(g);
    return g;

error:
    fprintf(stderr, "** file has invalid format **\n");
    free_graph(g);
    fclose(fp);
    free(line);
    return NULL;
}

struct graph *duplicate_graph(const struct graph *g)
{
    struct adjacency *adj1, *adj2;
    struct link *link;
    struct graph *dst;
    uint32_t i;

    if (!(dst = alloc_graph(g->num_nodes, g->flags)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
    {
        adj1 = &g->nodes[i];
        if (!adj1->links || !adj1->num_links) continue;

        if (!(link = xmalloc(sizeof(*link) * adj1->num_links)))
        {
            free_graph(dst);
            return NULL;
        }

        adj2 = &dst->nodes[i];
        adj2->ops       = adj1->ops;
        adj2->links     = link;
        adj2->num_links = adj1->num_links;
        adj2->max_links = adj1->num_links;
        adj2->hint      = adj1->hint;
        memcpy(link, adj1->links, sizeof(*link) * adj1->num_links);
    }

    return dst;
}

struct graph *transpose_graph(const struct graph *g)
{
    struct graph *dst;
    struct link *link;
    uint32_t i;

    if (!(dst = alloc_graph(g->num_nodes, g->flags & ~GRAPH_FLAGS_UNSORTED)))
        return NULL;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        _set_edge(dst, link->index, i, link->weight);
    }

    compress_graph_inline(dst);
    return dst;
}

struct graph *invert_graph(const struct graph *g, double max_weight, int self_loops)
{
    struct graph *dst;
    struct link *link;
    uint32_t i, j;

    if (!(dst = alloc_graph(g->num_nodes, g->flags & ~GRAPH_FLAGS_UNSORTED)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
    {
        GRAPH_FOR_EACH_LINK_ANY(g, &g->nodes[i], j, link)
        {
            float weight = max_weight;

            if (link)
            {
                weight -= link->weight;
                if (weight <= 0.0) continue;
            }

            if (i == j && !self_loops) continue;
            _set_edge(dst, i, j, weight);
        }
    }

    compress_graph_inline(dst);
    return dst;
}

struct graph *filter_graph_labels(const struct graph *g, uint32_t *labels)
{
    struct graph *dst;
    struct link *link;
    uint32_t i;

    if (!(dst = alloc_graph(g->num_nodes, g->flags & ~GRAPH_FLAGS_UNSORTED)))
        return NULL;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        if (labels[i] == labels[link->index]) continue;
        _add_edge(dst, i, link->index, link->weight);
    }

    compress_graph_inline(dst);
    return dst;
}

struct graph *filter_graph_weights(const struct graph *g, float min, float max)
{
    struct graph *dst;
    struct link *link;
    uint32_t i;

    if (!(dst = alloc_graph(g->num_nodes, g->flags & ~GRAPH_FLAGS_UNSORTED)))
        return NULL;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        if (link->weight < min) continue;
        if (link->weight > max) continue;
        _add_edge(dst, i, link->index, link->weight);
    }

    compress_graph_inline(dst);
    return dst;
}

struct graph *clamp_graph(const struct graph *g, float min, float max)
{
    struct graph *dst;
    struct link *link;
    uint32_t i;

    if (!(dst = alloc_graph(g->num_nodes, g->flags & ~GRAPH_FLAGS_UNSORTED)))
        return NULL;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        float weight = link->weight;
        if (weight < min) weight = min;
        else if (weight > max) weight = max;
        _add_edge(dst, i, link->index, weight);
    }

    compress_graph_inline(dst);
    return dst;
}

int clamp_graph_inplace(struct graph *g, float min, float max)
{
    struct link *link;
    uint32_t i;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        float weight = link->weight;
        if (weight < min) link->weight = min;
        else if (weight > max) link->weight = max;
    }

    return 1;
}

int resize_graph_inplace(struct graph *g, uint32_t num_nodes)
{
    const struct adjacency_ops *ops = &adjacency_ops_sorted;
    struct adjacency *nodes;
    uint32_t i;

    if (g->flags & GRAPH_FLAGS_UNSORTED)
        ops = &adjacency_ops_unsorted;

    if (num_nodes < g->num_nodes)
    {
        for (i = 0; i < num_nodes; i++)
        {
            struct adjacency *adj = sort_adjacency(&g->nodes[i]);
            while (adj->num_links > 0 && adj->links[adj->num_links - 1].index >= num_nodes)
                adj->num_links--;
        }

        for (; i < g->num_nodes; i++)
            free_adjacency(&g->nodes[i]);

        if (g->nodes != g->buffer)
        {
            if (num_nodes)
            {
                if ((nodes = xrealloc(g->nodes, sizeof(*nodes) * num_nodes)))
                    g->nodes = nodes;
            }
            else
            {
                free(g->nodes);
                g->nodes = NULL;
            }
        }

        g->num_nodes = num_nodes;
    }

    if (num_nodes > g->num_nodes)
    {
        if (g->nodes != g->buffer)
        {
            if (!(nodes = xrealloc(g->nodes, sizeof(*nodes) * num_nodes)))
                return 0;
        }
        else
        {
            if (!(nodes = xmalloc(sizeof(*nodes) * num_nodes))) return 0;
            memcpy(nodes, g->nodes, sizeof(*nodes) * g->num_nodes);
        }

        for (i = g->num_nodes; i < num_nodes; i++)
            init_adjacency(&nodes[i], ops);

        g->num_nodes = num_nodes;
        g->nodes = nodes;
    }

    return 1;
}

void assign_graph(struct graph *g, struct graph *src)
{
    struct adjacency *adj1, *adj2;
    uint32_t i;

    assert(g->num_nodes == src->num_nodes);
    /* Please note that userdata is not assigned. */

    g->flags = src->flags;
    for (i = 0; i < g->num_nodes; i++)
    {
        adj1 = &g->nodes[i];
        adj2 = &src->nodes[i];

        free_adjacency(adj1);
        *adj1 = *adj2;
        init_adjacency(adj2, &adjacency_ops_unsorted);
    }

    free_graph(src);
}

void print_graph(FILE *file, const struct graph *g)
{
    struct link *link;
    uint32_t i;

    fprintf(file, "%u\n", g->num_nodes);
    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        if (link->weight == 1.0) fprintf(file, "%u %u\n", i + 1, link->index + 1);
        else fprintf(file, "%u %u %f\n", i + 1, link->index + 1, link->weight);
    }
    fprintf(file, "\n");
}

void compress_graph_inline(struct graph *g)
{
    struct adjacency *adj;
    struct link *link;
    uint32_t i;

    for (i = 0; i < g->num_nodes; i++)
    {
        adj = sort_adjacency(&g->nodes[i]);
        if (adj->num_links >= adj->max_links) continue;
        if (adj->num_links)
        {
            if (!(link = xrealloc(adj->links, sizeof(*link) * adj->num_links)))
                continue;

            adj->max_links = adj->num_links;
            adj->links     = link;
        }
        else
        {
            free_adjacency(adj);
        }
    }
}

void free_graph(struct graph *g)
{
    uint32_t i;

    if (!g) return;
    for (i = 0; i < g->num_nodes; i++)
        free_adjacency(&g->nodes[i]);

    if (g->nodes != g->buffer)
        free(g->nodes);
    free(g);
}

int graph_has_edge(const struct graph *g, uint32_t start, uint32_t end)
{
    if (start >= g->num_nodes || end >= g->num_nodes) return 0;
    return _get_link(&g->nodes[start], end, 0) != NULL;
}

float graph_get_edge(const struct graph *g, uint32_t start, uint32_t end)
{
    struct link *link;
    if (start >= g->num_nodes || end >= g->num_nodes) return 0.0;
    if ((link = _get_link(&g->nodes[start], end, 0))) return link->weight;
    return 0.0;
}

uint64_t graph_get_edges(const struct graph *g, uint32_t *edges, float *weights, uint64_t max_edges)
{
    uint64_t count = 0;
    struct link *link;
    uint32_t i;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        if (!(g->flags & GRAPH_FLAGS_DIRECTED) && link->index < i) continue;
        if (count++ >= max_edges) continue;
        if (edges)
        {
            *edges++ = i;
            *edges++ = link->index;
        }
        if (weights)
        {
            *weights++ = link->weight;
        }
    }

    return count;
}

void graph_get_weights(const struct graph *g, float *weights)
{
    struct link *link;
    uint32_t i, j;
    float *val = weights;

    for (i = 0; i < g->num_nodes; i++)
    {
        GRAPH_FOR_EACH_LINK_ANY(g, &g->nodes[i], j, link)
        {
            *val++ = link ? link->weight : 0.0;
        }
    }
}

void graph_set_edge(struct graph *g, uint32_t start, uint32_t end, float weight)
{
    struct link *link;
    if (start >= g->num_nodes || end >= g->num_nodes) return;
    if ((link = _get_link(&g->nodes[start], end, 1))) link->weight = weight;
    if ((g->flags & GRAPH_FLAGS_DIRECTED) || start == end) return;
    if ((link = _get_link(&g->nodes[end], start, 1))) link->weight = weight;
}

void graph_set_edges(struct graph *g, uint32_t *edges, float *weights, uint64_t num_edges)
{
    while (num_edges--)
    {
        graph_set_edge(g, edges[0], edges[1], weights[0]);
        edges += 2;
        weights++;
    }
}

void graph_add_edge(struct graph *g, uint32_t start, uint32_t end, float weight)
{
    struct link *link;
    if (start >= g->num_nodes || end >= g->num_nodes) return;
    if ((link = _get_link(&g->nodes[start], end, 1))) link->weight += weight;
    if ((g->flags & GRAPH_FLAGS_DIRECTED) || start == end) return;
    if ((link = _get_link(&g->nodes[end], start, 1))) link->weight += weight;
}

void graph_add_edges(struct graph *g, uint32_t *edges, float *weights, uint64_t num_edges)
{
    while (num_edges--)
    {
        graph_add_edge(g, edges[0], edges[1], weights[0]);
        edges += 2;
        weights++;
    }
}

void graph_min_edge(struct graph *g, uint32_t start, uint32_t end, float weight)
{
    struct link *link;
    if (start >= g->num_nodes || end >= g->num_nodes) return;
    if ((link = _get_link(&g->nodes[start], end, 1))) link->weight = MIN(link->weight, weight);
    if ((g->flags & GRAPH_FLAGS_DIRECTED) || start == end) return;
    if ((link = _get_link(&g->nodes[end], start, 1))) link->weight = MIN(link->weight, weight);
}

void graph_min_edges(struct graph *g, uint32_t *edges, float *weights, uint64_t num_edges)
{
    while (num_edges--)
    {
        graph_min_edge(g, edges[0], edges[1], weights[0]);
        edges += 2;
        weights++;
    }
}

void graph_max_edge(struct graph *g, uint32_t start, uint32_t end, float weight)
{
    struct link *link;
    if (start >= g->num_nodes || end >= g->num_nodes) return;
    if ((link = _get_link(&g->nodes[start], end, 1))) link->weight = MAX(link->weight, weight);
    if ((g->flags & GRAPH_FLAGS_DIRECTED) || start == end) return;
    if ((link = _get_link(&g->nodes[end], start, 1))) link->weight = MAX(link->weight, weight);
}

void graph_max_edges(struct graph *g, uint32_t *edges, float *weights, uint64_t num_edges)
{
    while (num_edges--)
    {
        graph_max_edge(g, edges[0], edges[1], weights[0]);
        edges += 2;
        weights++;
    }
}

void graph_del_edge(struct graph *g, uint32_t start, uint32_t end)
{
    struct link *link;
    if (start >= g->num_nodes || end >= g->num_nodes) return;
    if ((link = _get_link(&g->nodes[start], end, 0))) _del_link(&g->nodes[start], link);
    if ((g->flags & GRAPH_FLAGS_DIRECTED) || start == end) return;
    if ((link = _get_link(&g->nodes[end], start, 0))) _del_link(&g->nodes[end], link);
}

void graph_del_edges(struct graph *g, uint32_t *edges, uint64_t num_edges)
{
    while (num_edges--)
    {
        graph_del_edge(g, edges[0], edges[1]);
        edges += 2;
    }
}

uint64_t graph_count_links(const struct graph *g)
{
    uint64_t count = 0;
    uint32_t i;

    for (i = 0; i < g->num_nodes; i++)
        count += g->nodes[i].num_links;

    /* for undirected graphs divide by two */
    if (!(g->flags & GRAPH_FLAGS_DIRECTED))
    {
        count += graph_count_links_diagonal(g);
        assert(!(count & 1));
        count /= 2;
    }

    return count;
}

/* similar to graph_count_links, but ignore g->flags */
uint64_t graph_count_links_directed(const struct graph *g)
{
    uint64_t count = 0;
    uint32_t i;

    for (i = 0; i < g->num_nodes; i++)
        count += g->nodes[i].num_links;

    return count;
}

uint64_t graph_count_links_diagonal(const struct graph *g)
{
    uint64_t count = 0;
    uint32_t i;

    for (i = 0; i < g->num_nodes; i++)
        if (_get_link(&g->nodes[i], i, 0)) count++;

    return count;
}

double graph_sum_weights(const struct graph *g)
{
    struct link *link;
    double sum = 0.0;
    uint32_t i;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        sum += link->weight;
    }

    /* for undirected graphs divide by two */
    if (!(g->flags & GRAPH_FLAGS_DIRECTED))
    {
        for (i = 0; i < g->num_nodes; i++)
        {
            if (!(link = _get_link(&g->nodes[i], i, 0))) continue;
            sum += link->weight;
        }

        sum /= 2.0;
    }

    return sum;
}

uint32_t graph_out_degree(const struct graph *g, uint32_t index)
{
    struct adjacency *adj;

    if (index >= g->num_nodes) return 0;

    adj = &g->nodes[index];
    return adj->num_links;
}

uint32_t *graph_out_degrees(const struct graph *g)
{
    uint32_t *degree;
    uint32_t i;

    if (!(degree = xmalloc(sizeof(*degree) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
        degree[i] = g->nodes[i].num_links;

    return degree;
}

uint32_t graph_out_edges(const struct graph *g, uint32_t index, uint32_t *edges, float *weights, uint32_t max_edges)
{
    uint32_t count = 0;
    struct link *link;

    if (index >= g->num_nodes) return 0;

    GRAPH_FOR_EACH_LINK(&g->nodes[index], link)
    {
        if (count++ >= max_edges) continue;
        if (edges) *edges++ = link->index;
        if (weights) *weights++ = link->weight;
    }

    return count;
}

uint32_t graph_in_degree(const struct graph *g, uint32_t index)
{
    uint32_t i, degree = 0;

    if (index >= g->num_nodes) return 0;

    for (i = 0; i < g->num_nodes; i++)
        if (_get_link(&g->nodes[i], index, 0)) degree++;

    return degree;
}

uint32_t *graph_in_degrees(const struct graph *g)
{
    struct link *link;
    uint32_t *degree;
    uint32_t i;

    if (!(degree = xmalloc(sizeof(*degree) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
        degree[i] = 0;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        degree[link->index]++;
    }

    return degree;
}

double *graph_out_weights(const struct graph *g)
{
    struct link *link;
    double *weight;
    uint32_t i;

    if (!(weight = xmalloc(sizeof(*weight) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
        weight[i] = 0.0;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        weight[i] += link->weight;
    }

    return weight;
}

double *graph_in_weights(const struct graph *g)
{
    struct link *link;
    double *weight;
    uint32_t i;

    if (!(weight = xmalloc(sizeof(*weight) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
        weight[i] = 0.0;

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        weight[link->index] += link->weight;
    }

    return weight;
}

double *graph_degree_anomalies(const struct graph *g)
{
    struct link *link;
    double *result;
    uint32_t i;

    assert(!(g->flags & GRAPH_FLAGS_DIRECTED));

    if (!(result = xmalloc(sizeof(*result) * g->num_nodes)))
        return NULL;

    for (i = 0; i < g->num_nodes; i++)
    {
        uint64_t value = 0;
        uint32_t local = 0;

        GRAPH_FOR_EACH_LINK(&g->nodes[i], link)
        {
            struct adjacency *adj = &g->nodes[link->index];
            value += adj->num_links;
            local++;
        }

        if (local) result[i] = local - (double)value / local;
        else result[i] = 0.0;
    }

    return result;
}

double *graph_weight_anomalies(const struct graph *g)
{
    struct link *link;
    double *weight;
    double *result;
    uint32_t i;

    assert(!(g->flags & GRAPH_FLAGS_DIRECTED));

    if (!(weight = graph_out_weights(g)))
        return NULL;

    if (!(result = xmalloc(sizeof(*result) * g->num_nodes)))
    {
        free(weight);
        return NULL;
    }

    for (i = 0; i < g->num_nodes; i++)
    {
        double value = 0;

        GRAPH_FOR_EACH_LINK(&g->nodes[i], link)
        {
            value += link->weight * weight[link->index];
        }

        if (weight[i] != 0.0) result[i] = weight[i] - value / weight[i];
        else result[i] = 0.0;
    }

    free(weight);
    return result;
}

struct adjacency *sort_adjacency(struct adjacency *adj)
{
    if (adj->ops != &adjacency_ops_sorted)
    {
        adj->ops = &adjacency_ops_sorted;
        adj->ops->init(adj);
    }
    return adj;
}
