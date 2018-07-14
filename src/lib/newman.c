/*
 * Minimal graph library
 * Modularity maximization/minimization functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#define _GNU_SOURCE
#include <inttypes.h>
#include "graph.h"
#include "internal.h"

/*
 * Implementation based on:
 *
 * - M. E. J. Newman. 2004. Fast algorithm for detecting community
 *   structure in networks. Phys. Rev. E 69, 6 (2004), 5.
 *   https://doi.org/10.1103/PhysRevE.69.066133
 *
 * - S. Lackner, A. Spitz, M. Weidemüller, and M. Gertz. 2018. Efficient
 *   Anti-community Detection in Complex Networks. Proceedings of the 30th
 *   International Conference on Scientific and Statistical Database
 *   Management (SSDBM '18), July 9-11, 2018, Bozen-Bolzano, Italy. ACM,
 *   NY, USA. https://doi.org/10.1145/3221269.3221289
 *
 */

struct newman_merge
{
    uint32_t index0, index1;
    double delta_mod;
};

struct newman_context
{
    uint32_t flags;
    int (*sort_func)(const void *, const void *, void *);

    uint32_t *labels;
    double norm;

    struct
    {
        double *weights;
        uint32_t *lookup;
    } out, in;
};

static uint32_t *alloc_lookup(struct newman_context *context, struct graph *g, double *weights)
{
    uint32_t *lookup = xmalloc(sizeof(*lookup) * g->num_nodes);
    uint32_t i;

    for (i = 0; i < g->num_nodes; i++) lookup[i] = i;                           /* O(n) */
    qsort_r(lookup, g->num_nodes, sizeof(lookup[0]),                            /* O(n log(n)) */
            context->sort_func, weights);

    return lookup;
}

static void resort_lookup(struct newman_context *context, struct graph *g, uint32_t *lookup, double *weights, uint32_t index)
{
    uint32_t i;

    if (!lookup) return;

    for (i = 0; i < g->num_nodes; i++)                                          /* O(n) */
        if (lookup[i] == index) break;
    assert(i < g->num_nodes);
    resort_r(lookup, g->num_nodes, sizeof(lookup[0]), i,
             context->sort_func, weights);
}

static inline double init_metric(struct newman_context *context)
{
    return (context->flags & NEWMAN_FLAGS_MINIMIZE) ? INFINITY : -INFINITY;
}

/* maximization: value1 > value2,
 * minimization: value1 < value2 */
static inline int compare_metric_gt(struct newman_context *context, double value1, double value2)
{
    return (context->flags & NEWMAN_FLAGS_MINIMIZE) ? (value1 < value2) : (value1 > value2);
}

/* maximization: value1 >= value2,
 * minimization: value1 <= value2 */
static inline int compare_metric_ge(struct newman_context *context, double value1, double value2)
{
    return (context->flags & NEWMAN_FLAGS_MINIMIZE) ? (value1 <= value2) : (value1 >= value2);
}

static inline float modularity_get_edge(struct newman_context *context, uint32_t start, uint32_t end)
{
    return - context->out.weights[start] * context->in.weights[end] / context->norm;
}

static int sort_index_by_weight_inc(const void *a, const void *b, void *userdata)
{
    const double *weights = userdata;
    const uint32_t *ia = a, *ib = b;
    return COMPARE(weights[*ia], weights[*ib]);
}

static int sort_index_by_weight_dec(const void *a, const void *b, void *userdata)
{
    const double *weights = userdata;
    const uint32_t *ia = a, *ib = b;
    return COMPARE(weights[*ib], weights[*ia]);
}

/* COMPLEXITY: O(m) */
static int get_merge_edge_directed(struct newman_context *context, struct graph *g, struct newman_merge *merge)
{
    struct link *link;
    float weight;
    uint32_t i;
    int found = 0;

    assert(g->flags & GRAPH_FLAGS_DIRECTED);
    merge->delta_mod = init_metric(context);

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        if (link->index == i) continue;

        weight  = link->weight + modularity_get_edge(context, i, link->index);
        weight += graph_get_edge(g, link->index, i) + modularity_get_edge(context, link->index, i);
        if (compare_metric_gt(context, weight, merge->delta_mod))
        {
            merge->index0 = i;
            merge->index1 = link->index;
            merge->delta_mod = weight;
            found = 1;
        }
    }

    return found;
}

/* COMPLEXITY: O(m) */
static int get_merge_edge_undirected(struct newman_context *context, struct graph *g, struct newman_merge *merge)
{
    struct link *link;
    float weight;
    uint32_t i;
    int found = 0;

    assert(!(g->flags & GRAPH_FLAGS_DIRECTED));
    merge->delta_mod = init_metric(context);

    GRAPH_FOR_EACH_EDGE_UNDIRECTED(g, i, link)
    {
        if (link->index == i) continue;

        weight = 2.0 * (link->weight + modularity_get_edge(context, i, link->index));
        if (compare_metric_gt(context, weight, merge->delta_mod))
        {
            merge->index0 = i;
            merge->index1 = link->index;
            merge->delta_mod = weight;
            found = 1;
        }
    }

    return found;
}

/* COMPLEXITY: O(m log(n)) */
static int get_merge_noedge_undirected(struct newman_context *context, struct graph *g, struct newman_merge *merge)
{
    static uint32_t *state, num_states;
    double value, lower, upper = init_metric(context);
    uint32_t i, j, k, l;
    int found = 0;

    assert(!(g->flags & GRAPH_FLAGS_DIRECTED));

    if (!context->out.lookup)
    {
        context->out.lookup = alloc_lookup(context, g, context->out.weights);                /* O(n log(n)) */
        context->in.lookup  = context->out.lookup;
    }

    if (g->num_nodes != num_states)
    {
        state = xrealloc(state, sizeof(*state) * g->num_nodes);
        num_states = g->num_nodes;
    }

#define IGNORE(x) ({uint32_t __i = context->out.lookup[x]; context->labels[__i] != __i;})
#define WEIGHT(x) context->out.weights[context->out.lookup[x]]

    for (i = 0; i < g->num_nodes; i++)
    {
        if (IGNORE(i)) state[i] = g->num_nodes;
        else state[i] = i + 1;
    }

    merge->delta_mod = init_metric(context);
    for (i = 0; i < g->num_nodes - 1; i++)
    {
        lower = WEIGHT(i + 1) * WEIGHT(i + 1);
        for (j = 0; j <= i; j++)
        {
            for (; state[j] < g->num_nodes; state[j]++)
            {
                value = WEIGHT(j) * WEIGHT(state[j]);
                assert(compare_metric_ge(context, value, upper));

                if (compare_metric_gt(context, value, lower)) break;
                if (found && compare_metric_ge(context, value, lower)) break;  /* not better than existing */
                if (IGNORE(state[j])) continue;

                k = context->out.lookup[j];
                l = context->out.lookup[state[j]];
                if (!graph_has_edge(g, k, l))
                {
                    merge->index0 = k;
                    merge->index1 = l;
                    merge->delta_mod = 2.0 * modularity_get_edge(context, k, l);
                    lower = value;
                    found = 1;
                    break;
                }
            }
        }

        if (found) return 1;
        upper = lower;
    }

#undef IGNORE
#undef WEIGHT

    for (i = 0; i < g->num_nodes; i++)
        assert(state[i] >= g->num_nodes);

    return 0;
}

/* COMPLEXITY: O(n^2 log(n)) */
static int get_merge_any_directed(struct newman_context *context, struct graph *g, struct newman_merge *merge)
{
    uint32_t i, j;
    float weight;
    int found = 0;

    assert(g->flags & GRAPH_FLAGS_DIRECTED);
    merge->delta_mod = init_metric(context);

    if (context->flags & NEWMAN_FLAGS_FAST)
        fprintf(stderr, "NEWMAN_FLAGS_FAST not supported here\n");

    for (i = 0; i < g->num_nodes; i++)
    {
        if (context->labels[i] != i) continue;
        for (j = 0; j < i; j++)
        {
            if (context->labels[j] != j) continue;

            weight  = graph_get_edge(g, i, j) + modularity_get_edge(context, i, j);
            weight += graph_get_edge(g, j, i) + modularity_get_edge(context, j, i);
            if (compare_metric_gt(context, weight, merge->delta_mod))
            {
                merge->index0 = i;
                merge->index1 = j;
                merge->delta_mod = weight;
                found = 1;
            }
        }
    }

    return found;
}

/* COMPLEXITY: O(m + n) */
static int get_merge_any_undirected(struct newman_context *context, struct graph *g, struct newman_merge *merge)
{
    struct link *link;
    float weight;
    uint32_t i;
    int found;

    assert(!(g->flags & GRAPH_FLAGS_DIRECTED));

    if (!(found = get_merge_noedge_undirected(context, g, merge)))
        merge->delta_mod = init_metric(context);

    GRAPH_FOR_EACH_EDGE_UNDIRECTED(g, i, link)
    {
        if (link->index == i) continue;

        weight = 2.0 * (link->weight + modularity_get_edge(context, i, link->index));
        if (compare_metric_gt(context, weight, merge->delta_mod))
        {
            merge->index0 = i;
            merge->index1 = link->index;
            merge->delta_mod = weight;
            found = 1;
        }
    }

    return found;
}

static int (*newman_funcs[2][2][2])(struct newman_context *, struct graph *, struct newman_merge *) =
{
    {{ get_merge_any_undirected,        /* undirected, maximize, slow */
       get_merge_edge_undirected },     /* undirected, maximize, fast */
     { get_merge_any_undirected,        /* undirected, minimize, slow */
       get_merge_noedge_undirected }},  /* undirected, minimize, fast */
    {{ get_merge_any_directed,          /* directed, maximize, slow */
       get_merge_edge_directed },       /* directed, maximize, fast */
     { get_merge_any_directed,          /* directed, minimize, slow */
       get_merge_any_directed }},       /* directed, minimize, fast */      /* FIXME: Not implemented */
};

uint32_t *newman(struct graph *g, uint32_t flags, uint32_t num_clusters)
{
    int (*merge_func)(struct newman_context *, struct graph *, struct newman_merge *);
    struct newman_context context;
    struct newman_merge merge;
    uint32_t *labels = NULL;
    uint32_t num_iter, i;
    double best, mod = 0.0;

    merge_func = newman_funcs[(g->flags & GRAPH_FLAGS_DIRECTED) != 0]
                             [(flags & NEWMAN_FLAGS_MINIMIZE) != 0]
                             [(flags & NEWMAN_FLAGS_FAST) != 0];

    /* setup callback context */
    context.flags       = flags;
    context.sort_func   = (flags & NEWMAN_FLAGS_MINIMIZE) ? sort_index_by_weight_dec : sort_index_by_weight_inc;
    context.labels      = alloc_labels(g->num_nodes);
    context.norm        = graph_sum_weights(g);                                 /* O(m) */

    context.out.weights = graph_out_weights(g);                                 /* O(m) */
    context.out.lookup  = NULL;

    if (g->flags & GRAPH_FLAGS_DIRECTED)
    {
        context.in.weights = graph_in_weights(g);                               /* O(m) */
        context.in.lookup  = NULL;
    }
    else
    {
        context.in.weights = context.out.weights;
        context.in.lookup  = NULL;
        context.norm *= 2.0;
    }

    if (flags & NEWMAN_FLAGS_ANTIMODULARITY)
    {
        context.norm = g->num_nodes;
        assign_graph(g, square_graph(g));                                       /* O(mn log(n)) */
        fprintf(stderr, "Squared graph has %" PRIu64 "/%" PRIu64 " edges\n",
                graph_count_links_directed(g), (uint64_t)g->num_nodes * g->num_nodes);
    }

    /* initialize modularity */
    for (i = 0; i < g->num_nodes; i++)                                          /* O(n log(n)) */
        mod += graph_get_edge(g, i, i) + modularity_get_edge(&context, i, i);

    /* track best result */
    if (!num_clusters || num_clusters == count_labels(context.labels, g->num_nodes))
    {
        if (!labels) labels = alloc_labels(g->num_nodes);
        memcpy(labels, context.labels, sizeof(*labels) * g->num_nodes);
        best = mod;
    }
    else
    {
        best = init_metric(&context);
    }

    for (num_iter = 0; merge_func(&context, g, &merge); num_iter++)             /* O(MERGE) */
    {
        progress("Processing %u / %u", num_iter, g->num_nodes);

        /* update modularity and labels */
        assert(merge.index0 != merge.index1);
        assert(context.labels[merge.index0] == merge.index0);
        assert(context.labels[merge.index1] == merge.index1);

        mod += merge.delta_mod;
        for (i = 0; i < g->num_nodes; i++)                                      /* O(n) */
        {
            if (context.labels[i] == merge.index1)
                context.labels[i] = merge.index0;
        }

        /* update modularity matrix */
        graph_merge_nodes_add(g, merge.index0, merge.index1);                   /* O(n log(n)) + O(n^2) */

        /* update weights */
        context.out.weights[merge.index0] += context.out.weights[merge.index1];
        resort_lookup(&context, g, context.out.lookup, context.out.weights, merge.index0);    /* O(n) */
        context.out.weights[merge.index1] = 0.0;
        resort_lookup(&context, g, context.out.lookup, context.out.weights, merge.index1);    /* O(n) */

        if (g->flags & GRAPH_FLAGS_DIRECTED)
        {
            context.in.weights[merge.index0] += context.in.weights[merge.index1];
            resort_lookup(&context, g, context.in.lookup, context.in.weights, merge.index0);  /* O(n) */
            context.in.weights[merge.index1] = 0.0;
            resort_lookup(&context, g, context.in.lookup, context.in.weights, merge.index1);  /* O(n) */
        }

        /* track best result */
        if (compare_metric_gt(&context, mod, best) &&
            (!num_clusters || num_clusters == count_labels(context.labels, g->num_nodes)))
        {
            if (!labels) labels = alloc_labels(g->num_nodes);
            memcpy(labels, context.labels, sizeof(*labels) * g->num_nodes);
            best = mod;
        }

        /* in verbose mode print what we are doing */
        if (flags & NEWMAN_FLAGS_VERBOSE)
        {
            fprintf(stderr, "Merging community %u and %u (delta_mod = %f)\n",
                    merge.index0, merge.index1, merge.delta_mod);
            fprintf(stderr, "M([");
            for (i = 0; i < g->num_nodes; i++)
            {
                if (i) fprintf(stderr, ", ");
                fprintf(stderr, "%u", context.labels[i]);
            }
            fprintf(stderr, "]) = %f (unnormalized)\n", mod);
        }

        /* in strict mode, validate that the modularity is computed correctly */
        if (flags & NEWMAN_FLAGS_STRICT)
        {
            double new_mod = 0.0;
            for (i = 0; i < g->num_nodes; i++)
            {
                if (context.labels[i] != i) continue;
                new_mod += graph_get_edge(g, i, i) + modularity_get_edge(&context, i, i);
            }
            if (fabs(mod - new_mod) >= 1e-5)
            {
                fprintf(stderr, "Modularity mismatch: Expected %f, got %f\n", new_mod, mod);
                assert(0);
            }
            if (context.out.lookup)
            {
                for (i = 1; i < g->num_nodes; i++)
                {
                    assert(compare_metric_ge(&context, context.out.weights[context.out.lookup[i]],
                                                       context.out.weights[context.out.lookup[i - 1]]));
                    assert(compare_metric_ge(&context, context.out.weights[context.in.lookup[i]],
                                                       context.out.weights[context.in.lookup[i - 1]]));
                }
            }
        }
    }

    fprintf(stderr, "Modularity (unnormalized): %f\n", best);

    free(context.out.weights);
    free(context.out.lookup);
    if (g->flags & GRAPH_FLAGS_DIRECTED)
    {
        free(context.in.weights);
        free(context.in.lookup);
    }
    free_labels(context.labels);
    return labels;
}
