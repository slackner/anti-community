/*
 * Compute metrics for large graphs
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include <inttypes.h>
#include "graph.h"

#define METRIC_FUNC(name) \
    static double metric_ ## name(const uint32_t *labels_true, const uint32_t *labels_pred, const struct graph *g)

#define METRIC_WRAP(name) \
    METRIC_FUNC(name) { return labels_true ? name(labels_true, labels_pred, g->num_nodes) : NAN; }

METRIC_FUNC(nodes)          { return g->num_nodes; }
METRIC_FUNC(edges)          { return graph_count_links(g); }
METRIC_FUNC(diagonal)       { return graph_count_links_diagonal(g); }
METRIC_FUNC(communities)    { return count_labels(labels_pred, g->num_nodes); }
METRIC_FUNC(modularity)     { return modularity(g, labels_pred); }
METRIC_FUNC(antimodularity) { return antimodularity(g, labels_pred); }
METRIC_FUNC(density)
{
    double max_links;

    if (graph_count_links_diagonal(g))
    {
        fprintf(stderr, "FIXME: Density not implemented for graphs with selfloops\n");
        return NAN;
    }

    max_links = g->num_nodes * (double)(g->num_nodes - 1);
    if (!(g->flags & GRAPH_FLAGS_DIRECTED)) max_links /= 2.0;

    return graph_count_links(g) / max_links;
}

METRIC_WRAP(precision);
METRIC_WRAP(recall);
METRIC_WRAP(rand_index);
METRIC_WRAP(fowlkes_mallows);
METRIC_WRAP(jaccard);
METRIC_WRAP(f1_measure);
METRIC_WRAP(adj_rand_index);
METRIC_WRAP(norm_mutual_info);
METRIC_FUNC(adj_rand_index_comp)
{
    if (!labels_true) return NAN;
    return adj_rand_index_comp(labels_true, labels_pred, g);
}

#undef METRIC_FUNC
#undef METRIC_WRAP

struct metric_function
{
    const char *name;
    double (*func)(const uint32_t *, const uint32_t *, const struct graph *);
};

static const struct metric_function metrics[] =
{
    { "#nodes",                 metric_nodes                },
    { "#edges",                 metric_edges                },
    { "#diagonal",              metric_diagonal             },
    { "#communities",           metric_communities          },
    { "modularity",             metric_modularity           },
    { "antimodularity",         metric_antimodularity       },
    { "density",                metric_density              },
    { "precision",              metric_precision            },
    { "recall",                 metric_recall               },
    { "rand-index",             metric_rand_index           },
    { "fowlkes-mallows",        metric_fowlkes_mallows      },
    { "jaccard",                metric_jaccard              },
    { "f1-measure",             metric_f1_measure           },
    { "adj-rand-index",         metric_adj_rand_index       },
    { "adj-rand-index-comp",    metric_adj_rand_index_comp  },
    { "norm-mutual-info",       metric_norm_mutual_info     },
};

static void usage(const char *argv0)
{
    const char *ptr = strrchr(argv0, '/');
    if (ptr) argv0 = ptr + 1;

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: ./%s [--clusters=filename] [--groundtruth=filename] [graphfile]\n", argv0);
    fprintf(stderr, "\n");
    fprintf(stderr, "Compute modularity and anti-modularity for the graph\n");
    fprintf(stderr, "specified by the given file. When a ground-truth label\n");
    fprintf(stderr, "file is given, also compute cluster quality metrics.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Configuration:\n");
    fprintf(stderr, "  --clusters=filename  Specify the path to the clusters file\n");
    fprintf(stderr, "  --groundtruth=file   Specify the path to the ground-truth file\n");
    fprintf(stderr, "  --directed           Interpret as directed graph\n");
    fprintf(stderr, "  --invert             Invert the graph\n");
    fprintf(stderr, "  --unweighted         Ignore weights of the original graph\n");
    fprintf(stderr, "  --help               Display this help and exit\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    const char *filename = NULL;
    const char *clusters = NULL;
    const char *groundtruth = NULL;
    uint32_t *labels_pred = NULL;
    uint32_t *labels_true = NULL;
    struct graph *g;
    uint32_t flags = 0;
    int unweighted = 0;
    int invert = 0;
    double value;
    int ret = 1;
    int i;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--help")) usage(argv[0]);
        else if (!strcmp(argv[i], "--directed")) flags |= GRAPH_FLAGS_DIRECTED;
        else if (!strcmp(argv[i], "--invert")) invert = 1;
        else if (!strcmp(argv[i], "--unweighted")) unweighted = 1;
        else if (!strncmp(argv[i], "--clusters=", 11)) clusters = argv[i] + 11;
        else if (!strncmp(argv[i], "--groundtruth=", 14)) groundtruth = argv[i] + 14;
        else if (!filename) filename = argv[i];
        else fprintf(stderr, "Ignoring command line argument '%s'\n", argv[i]);
    }

    if (!filename || !clusters)
    {
        usage(argv[0]);
        return 1;
    }

    if (!(g = load_graph(filename, flags)))
    {
        fprintf(stderr, "Failed to load '%s'\n", filename);
        return 1;
    }

    if (groundtruth && !(labels_true = load_labels(groundtruth, g->num_nodes)))
    {
        fprintf(stderr, "Failed to load '%s'\n", groundtruth);
        goto error;
    }

    if (!(labels_pred = load_labels(clusters, g->num_nodes)))
    {
        fprintf(stderr, "Failed to load '%s'\n", clusters);
        goto error;
    }

    if (unweighted)
    {
        assign_graph(g, clamp_graph(g, 1.0, 1.0));
        fprintf(stderr, "Input has been converted to unweighted graph\n");
    }

    if (invert)
    {
        assign_graph(g, invert_graph(g, 1.0, 1));
        fprintf(stderr, "Inverted graph has %" PRIu64 "/%" PRIu64 " edges\n",
                graph_count_links_directed(g), (uint64_t)g->num_nodes * g->num_nodes);
    }

    for (i = 0; i < sizeof(metrics) / sizeof(metrics[0]); i++)
    {
        value = metrics[i].func(labels_true, labels_pred, g);
        if (!isnan(value)) printf("%s: %f\n", metrics[i].name, value);
    }

    ret = 0;

error:
    free_graph(g);
    free_labels(labels_true);
    free_labels(labels_pred);
    return ret;
}
