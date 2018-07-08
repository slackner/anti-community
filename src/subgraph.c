/*
 * Extract subgraph
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include <inttypes.h>
#include "graph.h"

static void print_comments(const char *filename)
{
    ssize_t read;
    size_t len = 0;
    char *line = NULL;
    FILE *fp;

    if (!(fp = fopen(filename, "r")))
        return;

    while ((read = getline(&line, &len, fp)) > 0)
    {
        if (line[read - 1] == '\n') line[read - 1] = 0;
        if (line[0] == '#' || line[0] == ';')
            printf("%s\n", line);
    }

    fclose(fp);
    free(line);
}

static void usage(const char *argv0)
{
    const char *ptr = strrchr(argv0, '/');
    if (ptr) argv0 = ptr + 1;

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: ./%s [--clusters=filename] [graphfile]\n", argv0);
    fprintf(stderr, "\n");
    fprintf(stderr, "Extract subgraph specified by the given files.\n");
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
    uint32_t *labels = NULL;
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

    if (clusters && !(labels = load_labels(clusters, g->num_nodes)))
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

    if (labels)
    {
        assign_graph(g, filter_graph_labels(g, labels));
        fprintf(stderr, "Filtered graph has %" PRIu64 "/%" PRIu64 " edges\n",
                graph_count_links_directed(g), (uint64_t)g->num_nodes * g->num_nodes);
    }

    print_comments(filename);
    print_graph(stdout, g);

    free_graph(g);
    free_labels(labels);
    return 0;
}
