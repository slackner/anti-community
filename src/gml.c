/*
 * Convert graph to GML file
 *
 * Copyright (c) 2018 Sebastian Lackner
 */

#include <inttypes.h>
#include "graph.h"

struct attr
{
    struct attr *next;
    const char  *filename;
    uint32_t    *labels;
    char         name[0];
};

static void print_graph_gml(FILE *file, const struct graph *g, struct attr *attrs)
{
    struct attr *a;
    struct link *link;
    uint32_t i;

    fprintf(file, "graph\n");
    fprintf(file, "[\n");

    if (g->flags & GRAPH_FLAGS_DIRECTED)
        fprintf(file, "  directed 1\n");

    for (i = 0; i < g->num_nodes; i++)
    {
        fprintf(file, "  node\n");
        fprintf(file, "  [\n");
        fprintf(file, "    id %u\n", i);
        for (a = attrs; a; a = a->next)
            fprintf(file, "    %s %u\n", a->name, a->labels[i]);
        fprintf(file, "  ]\n");
    }

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        if (!(g->flags & GRAPH_FLAGS_DIRECTED) && link->index < i) continue;
        fprintf(file, "  edge\n");
        fprintf(file, "  [\n");
        fprintf(file, "    source %u\n", i);
        fprintf(file, "    target %u\n", link->index);
        fprintf(file, "    value %f\n", link->weight);
        fprintf(file, "  ]\n");
    }
    fprintf(file, "]\n");
}

static void usage(const char *argv0)
{
    const char *ptr = strrchr(argv0, '/');
    if (ptr) argv0 = ptr + 1;

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: ./%s [--attr=filename] [graphfile]\n", argv0);
    fprintf(stderr, "\n");
    fprintf(stderr, "Convert graph to GML file format.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Configuration:\n");
    fprintf(stderr, "  --attr=filename      Specify the path to a node attribute file\n");
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
    struct attr *attrs = NULL, *a;
    struct graph *g = NULL;
    uint32_t flags = 0;
    int unweighted = 0;
    int invert = 0;
    const char *p;
    int ret = 1;
    int i;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--help")) usage(argv[0]);
        else if (!strcmp(argv[i], "--directed")) flags |= GRAPH_FLAGS_DIRECTED;
        else if (!strcmp(argv[i], "--invert")) invert = 1;
        else if (!strcmp(argv[i], "--unweighted")) unweighted = 1;
        else if (!strncmp(argv[i], "--", 2) && (p = strchr(argv[i] + 2, '=')))
        {
            if (!(a = malloc(offsetof(struct attr, name[p - (argv[i] + 2) + 1]))))
                goto error;

            a->next     = attrs;
            a->filename = p + 1;
            a->labels   = NULL;
            memcpy(a->name, argv[i] + 2, p - (argv[i] + 2));
            a->name[p - (argv[i] + 2)] = 0;

            attrs = a;
        }
        else if (!filename) filename = argv[i];
        else fprintf(stderr, "Ignoring command line argument '%s'\n", argv[i]);
    }

    if (!filename)
    {
        usage(argv[0]);
        goto error;
    }

    if (!(g = load_graph(filename, flags)))
    {
        fprintf(stderr, "Failed to load '%s'\n", filename);
        goto error;
    }

    for (a = attrs; a; a = a->next)
    {
        if (!(a->labels = load_labels(a->filename, g->num_nodes)))
        {
            fprintf(stderr, "Failed to load '%s'\n", a->filename);
            goto error;
        }
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

    print_graph_gml(stdout, g, attrs);
    ret = 0;

error:
    while ((a = attrs))
    {
        attrs = a->next;
        free_labels(a->labels);
        free(a);
    }
    free_graph(g);
    return ret;
}
