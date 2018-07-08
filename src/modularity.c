/*
 * Compute modularity for large graphs
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include <inttypes.h>
#include "graph.h"

static double node_count(const struct graph *g, const uint32_t *labels)
{
    return g->num_nodes;
}

static double edge_count(const struct graph *g, const uint32_t *labels)
{
    return graph_count_links(g);
}

static double diagonal_count(const struct graph *g, const uint32_t *labels)
{
    return graph_count_links_diagonal(g);
}

static double comm_count(const struct graph *g, const uint32_t *labels)
{
    return count_labels(labels, g->num_nodes);
}

static double density(const struct graph *g, const uint32_t *labels)
{
    double max_links;

    if (graph_count_links_diagonal(g))
    {
        fprintf(stderr, "FIXME: Not implemented for graphs with selfloops\n");
        return NAN;
    }

    max_links = g->num_nodes * (double)(g->num_nodes - 1);
    if (!(g->flags & GRAPH_FLAGS_DIRECTED)) max_links /= 2.0;

    return graph_count_links(g) / max_links;
}

struct metric_function
{
    const char *name;
    double (*func)(const struct graph *, const uint32_t *);
};

static const struct metric_function metrics[] =
{
    { "#nodes",         node_count },
    { "#edges",         edge_count },
    { "#diagonal",      diagonal_count },
    { "#communities",   comm_count },
    { "modularity",     modularity },
    { "antimodularity", antimodularity },
    { "density",        density },
};

static void usage(const char *argv0)
{
    const char *ptr = strrchr(argv0, '/');
    if (ptr) argv0 = ptr + 1;

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: ./%s [--clusters=filename] [graphfile]\n", argv0);
    fprintf(stderr, "\n");
    fprintf(stderr, "Compute modularity and anti-modularity for the graph\n");
    fprintf(stderr, "specified by the given file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Configuration:\n");
    fprintf(stderr, "  --clusters=filename  Specify the path to the clusters file\n");
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
    uint32_t *labels;
    struct graph *g;
    uint32_t flags = 0;
    int unweighted = 0;
    int invert = 0;
    int i;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--help")) usage(argv[0]);
        else if (!strcmp(argv[i], "--directed")) flags |= GRAPH_FLAGS_DIRECTED;
        else if (!strcmp(argv[i], "--invert")) invert = 1;
        else if (!strcmp(argv[i], "--unweighted")) unweighted = 1;
        else if (!strncmp(argv[i], "--clusters=", 11)) clusters = argv[i] + 11;
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

    if (!(labels = load_labels(clusters, g->num_nodes)))
    {
        fprintf(stderr, "Failed to load '%s'\n", clusters);
        free_graph(g);
        return 1;
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
        printf("%s: ", metrics[i].name);
        printf("%f\n", metrics[i].func(g, labels));
    }

    free_graph(g);
    free_labels(labels);
    return 0;
}
