/*
 * Convert graph to GML file
 *
 * Copyright (c) 2018 Sebastian Lackner
 */

#include <inttypes.h>
#include "graph.h"

static void usage(const char *argv0)
{
    const char *ptr = strrchr(argv0, '/');
    if (ptr) argv0 = ptr + 1;

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: ./%s [graphfile]\n", argv0);
    fprintf(stderr, "\n");
    fprintf(stderr, "Convert graph to GML file format.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Configuration:\n");
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
        else if (!filename) filename = argv[i];
        else fprintf(stderr, "Ignoring command line argument '%s'\n", argv[i]);
    }

    if (!filename)
    {
        usage(argv[0]);
        return 1;
    }

    if (!(g = load_graph(filename, flags)))
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

    print_graph_gml(stdout, g);
    free_graph(g);
    return 0;
}
