/*
 * Agglomerative hierarchical clustering algorithm
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include <inttypes.h>
#include "graph.h"

#ifndef NEWMAN_FLAGS
#define NEWMAN_FLAGS 0
#endif

static void usage(const char *argv0)
{
    const char *ptr = strrchr(argv0, '/');
    if (ptr) argv0 = ptr + 1;

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: ./%s [--verbose] [--strict] [graph]\n", argv0);
    fprintf(stderr, "\n");
    fprintf(stderr, "Agglomerative hierarchical clustering algorithm\n");
#if NEWMAN_FLAGS == NEWMAN_FLAGS_ANTIMODULARITY
    fprintf(stderr, "for anti-community detection / anti-modularity maximization.\n");
#elif NEWMAN_FLAGS == NEWMAN_FLAGS_MINIMIZE
    fprintf(stderr, "for anti-community detection / modularity minimization.\n");
#elif NEWMAN_FLAGS == 0
    fprintf(stderr, "for community detection / modularity maximization.\n");
#endif
    fprintf(stderr, "\n");
    fprintf(stderr, "Configuration:\n");
    fprintf(stderr, "  --clusters=N         Force a fixed number of clusters\n");
    fprintf(stderr, "  --directed           Interpret as directed graph\n");
    fprintf(stderr, "  --fast               Do not check non-existing edges\n");
    fprintf(stderr, "  --invert             Invert the graph\n");
    fprintf(stderr, "  --unweighted         Ignore weights of the original graph\n");
    fprintf(stderr, "  --verbose            Enable verbose output\n");
    fprintf(stderr, "  --strict             Enable self-checks (for debugging)\n");
    fprintf(stderr, "  --help               Display this help and exit\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    const char *filename = NULL;
    uint32_t num_clusters = 0;
    uint32_t *labels;
    struct graph *g;
    int newman_flags = NEWMAN_FLAGS;
    int graph_flags = 0;
    int unweighted = 0;
    int invert = 0;
    int i;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--help")) usage(argv[0]);
        else if (!strcmp(argv[i], "--verbose")) newman_flags |= NEWMAN_FLAGS_VERBOSE;
        else if (!strcmp(argv[i], "--strict")) newman_flags |= NEWMAN_FLAGS_STRICT;
        else if (!strcmp(argv[i], "--directed")) graph_flags |= GRAPH_FLAGS_DIRECTED;
        else if (!strcmp(argv[i], "--fast")) newman_flags |= NEWMAN_FLAGS_FAST;
        else if (!strcmp(argv[i], "--invert")) invert = 1;
        else if (!strcmp(argv[i], "--unweighted")) unweighted = 1;
        else if (!strncmp(argv[i], "--clusters=", 11)) num_clusters = atoi(argv[i] + 11);
        else if (!filename) filename = argv[i];
        else fprintf(stderr, "Ignoring command line argument '%s'\n", argv[i]);
    }

    if (!filename)
    {
        usage(argv[0]);
        return 1;
    }

    if (!(g = load_graph(filename, graph_flags)))
    {
        fprintf(stderr, "Failed to load '%s'\n", filename);
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

    if ((labels = newman(g, newman_flags, num_clusters)))
    {
        /* We cannot compute the modularity here, G has been modified. */
        print_labels(stdout, labels, g->num_nodes);
        free_labels(labels);
    }

    free_graph(g);
    return 0;
}
