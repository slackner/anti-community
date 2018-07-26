/*
 * Minimal graph library
 * Label vector manipulation functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#define _GNU_SOURCE
#include <math.h>
#include "graph.h"
#include "internal.h"

uint32_t *alloc_labels(uint32_t num_nodes)
{
    uint32_t *labels;
    uint32_t i;

    if (!(labels = xmalloc(sizeof(*labels) * num_nodes)))
        return NULL;

    for (i = 0; i < num_nodes; i++)
        labels[i] = i;

    return labels;
}

uint32_t *decode_labels(const uint32_t *indices, uint32_t num_nodes)
{
    uint32_t *labels;
    uint32_t label, i, j;

    if (!(labels = xmalloc(sizeof(*labels) * num_nodes)))
        return NULL;

    for (i = 0; i < num_nodes; i++)
        labels[i] = ~0U;

    for (i = 0; i < num_nodes; i++)
    {
        if (labels[i] != ~0U) continue;
        for (j = i; labels[j] == ~0U; j = indices[j])
        {
            labels[j] = i;
            if (indices[j] >= num_nodes) break;
        }

        if ((label = labels[j]) == i) continue;
        for (j = i; labels[j] == i; j = indices[j])
        {
            labels[j] = label;
            if (indices[j] >= num_nodes) break;
        }
    }

    return labels;
}

uint32_t *load_labels(const char *filename, uint32_t num_nodes)
{
    uint32_t label, pos = 0;
    char *ptr, *line = NULL;
    uint32_t *labels;
    ssize_t read;
    size_t len = 0;
    int offset;
    FILE *fp;

    if (!(labels = xmalloc(sizeof(*labels) * num_nodes)))
        return NULL;

    if (!(fp = fopen(filename, "r")))
    {
        free(labels);
        return NULL;
    }

    while ((read = getline(&line, &len, fp)) > 0)
    {
        if (line[read - 1] == '\n') line[read - 1] = 0;
        if (!line[0] || line[0] == '#' || line[0] == ';') continue;
        ptr = line;
        while (*ptr)
        {
            if (sscanf(ptr, "%u%n", &label, &offset) < 1) goto error;
            if (pos >= num_nodes) goto error;
            labels[pos++] = label;
            ptr += offset;
        }
    }

    fclose(fp);
    free(line);
    return labels;

error:
    fprintf(stderr, "** file has invalid format **\n");
    free(labels);
    fclose(fp);
    free(line);
    return NULL;
}

void print_labels(FILE *file, uint32_t *labels, uint32_t num_nodes)
{
    uint32_t i;

    /* FIXME: add some linebreaks */
    fprintf(file, "# Labels\n");
    for (i = 0; i < num_nodes; i++)
    {
        if (i) fprintf(file, " ");
        fprintf(file, "%u", labels[i]);
    }
    fprintf(file, "\n");
}

uint32_t count_labels(const uint32_t *labels, uint32_t num_nodes)
{
    uint32_t *lookup;
    uint32_t i, count = 0;

    if (!(lookup = xmalloc(sizeof(*lookup) * num_nodes)))
        return ~0U;

    for (i = 0; i < num_nodes; i++) lookup[i] = i;
    qsort_r(lookup, num_nodes, sizeof(lookup[0]), _sort_index_by_label, (void *)labels);

    for (i = 0; i < num_nodes;)
    {
        uint32_t label = labels[lookup[i]];
        while (i < num_nodes && labels[lookup[i]] == label) i++;
        count++;
    }

    free(lookup);
    return count;
}

uint32_t *simplify_labels(const uint32_t *labels, uint32_t num_nodes, uint32_t *count_out)
{
    uint32_t *output, *lookup;
    uint32_t i, count = 0;

    if (!(output = xmalloc(sizeof(*output) * num_nodes)))
        return NULL;

    if (!(lookup = xmalloc(sizeof(*lookup) * num_nodes)))
    {
        free(output);
        return NULL;
    }

    for (i = 0; i < num_nodes; i++) lookup[i] = i;
    qsort_r(lookup, num_nodes, sizeof(lookup[0]), _sort_index_by_label, (void *)labels);

    for (i = 0; i < num_nodes;)
    {
        uint32_t label = labels[lookup[i]];
        while (i < num_nodes && labels[lookup[i]] == label)
            output[lookup[i++]] = count;
        count++;
    }

    if (count_out)
        *count_out = count;

    free(lookup);
    return output;
}

uint32_t simplify_labels_inplace(uint32_t *labels, uint32_t num_nodes)
{
    uint32_t *lookup;
    uint32_t i, count = 0;

    if (!(lookup = xmalloc(sizeof(*lookup) * num_nodes)))
        return ~0U;

    for (i = 0; i < num_nodes; i++) lookup[i] = i;
    qsort_r(lookup, num_nodes, sizeof(lookup[0]), _sort_index_by_label, (void *)labels);

    for (i = 0; i < num_nodes;)
    {
        uint32_t label = labels[lookup[i]];
        while (i < num_nodes && labels[lookup[i]] == label)
            labels[lookup[i++]] = count;
        count++;
    }

    free(lookup);
    return count;
}

void split_labels(uint32_t *labels1, const uint32_t *labels2, uint32_t num_nodes)
{
    struct sort_index_by_label2_context context;
    uint32_t *lookup;
    uint32_t i, count = 0;

    if (!(lookup = xmalloc(sizeof(*lookup) * num_nodes)))
        return;

    context.labels1 = labels1;
    context.labels2 = labels2;
    for (i = 0; i < num_nodes; i++) lookup[i] = i;
    qsort_r(lookup, num_nodes, sizeof(lookup[0]), _sort_index_by_label2, &context);

    for (i = 0; i < num_nodes;)
    {
        uint32_t index = lookup[i];
        uint32_t label1 = labels1[index];
        uint32_t label2 = labels2[index];
        while (i < num_nodes)
        {
            index = lookup[i];
            if (labels1[index] != label1 || labels2[index] != label2) break;
            labels1[index] = count;
            i++;
        }
        count++;
    }

    free(lookup);
}

struct graph *intersection_matrix(const uint32_t *labels1, const uint32_t *labels2, uint32_t num_nodes)
{
    uint32_t *labels1_simple, *labels2_simple;
    uint32_t labels1_count, labels2_count;
    uint32_t max_labels;
    struct graph *g;
    uint32_t i;

    if (!(labels1_simple = simplify_labels(labels1, num_nodes, &labels1_count)))
        return NULL;

    if (!(labels2_simple = simplify_labels(labels2, num_nodes, &labels2_count)))
    {
        free(labels1_simple);
        return NULL;
    }

    max_labels = MAX(labels1_count, labels2_count);
    if ((g = alloc_graph(max_labels, GRAPH_FLAGS_DIRECTED)))
    {
        for (i = 0; i < num_nodes; i++)
            _add_edge(g, labels1_simple[i], labels2_simple[i], 1.0);
    }

    free(labels1_simple);
    free(labels2_simple);
    return g;
}

/* a: true positive, b: false positive, c: false negative, d: true negative */
int confusion_matrix(double *a_out, double *b_out, double *c_out, double *d_out,
                     const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes)
{
    double a = 0.0, b = 0.0, c = 0.0;
    struct link *link;
    double *weights;
    struct graph *g;
    uint32_t i;

    if (!(g = intersection_matrix(labels_true, labels_pred, num_nodes)))
        return 0;

    GRAPH_FOR_EACH_EDGE(g, i, link)
        a += link->weight * (link->weight - 1);

    if (!(weights = graph_in_weights(g)))
    {
        free_graph(g);
        return 0;
    }

    for (i = 0; i < g->num_nodes; i++)
        b += weights[i] * (weights[i] - 1);

    free(weights);
    if (!(weights = graph_out_weights(g)))
    {
        free_graph(g);
        return 0;
    }

    for (i = 0; i < g->num_nodes; i++)
        c += weights[i] * (weights[i] - 1);

    free(weights);

    b -= a;
    c -= a;

    if (a_out) *a_out = a / 2.0;
    if (b_out) *b_out = b / 2.0;
    if (c_out) *c_out = c / 2.0;
    if (d_out) *d_out = ((double)num_nodes * (num_nodes - 1) - a - b - c) / 2.0;

    free_graph(g);
    return 1;
}

double precision(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes)
{
    double a, b, c;

    if (!confusion_matrix(&a, &b, &c, NULL, labels_true, labels_pred, num_nodes)) return NAN;
    if (b == 0.0 && c == 0.0) return 1.0;
    if (a == 0.0) return 0.0;  /* FIXME */

    return a / (a + b);
}

double recall(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes)
{
    double a, b, c;

    if (!confusion_matrix(&a, &b, &c, NULL, labels_true, labels_pred, num_nodes)) return NAN;
    if (b == 0.0 && c == 0.0) return 1.0;
    if (a == 0.0) return 0.0;  /* FIXME */

    return a / (a + c);
}

double rand_index(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes)
{
    double a, b, c, d;

    if (!confusion_matrix(&a, &b, &c, &d, labels_true, labels_pred, num_nodes)) return NAN;
    if (b == 0.0 && c == 0.0) return 1.0;

    return (a + d) / (a + b + c + d);
}

double fowlkes_mallows(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes)
{
    double a, b, c;

    if (!confusion_matrix(&a, &b, &c, NULL, labels_true, labels_pred, num_nodes)) return NAN;
    if (b == 0.0 && c == 0.0) return 1.0;
    if (a == 0.0) return 0.0;  /* FIXME */

    return a / sqrt((a + b) * (a + c));
}

double jaccard(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes)
{
    double a, b, c;

    if (!confusion_matrix(&a, &b, &c, NULL, labels_true, labels_pred, num_nodes)) return NAN;
    if (b == 0.0 && c == 0.0) return 1.0;

    return a / (a + b + c);
}

double f1_measure(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes)
{
    double a, b, c;

    if (!confusion_matrix(&a, &b, &c, NULL, labels_true, labels_pred, num_nodes)) return NAN;
    if (b == 0.0 && c == 0.0) return 1.0;

    return (2.0 * a) / (2.0 * a + b + c);
}

double adjusted_rand_index(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes)
{
    double a, b, c, d;
    double n2, tmp;

    if (!confusion_matrix(&a, &b, &c, &d, labels_true, labels_pred, num_nodes)) return NAN;
    if (b == 0.0 && c == 0.0) return 1.0;

    n2  = a + b + c + d;
    tmp = (a + b) * (a + c) + (c + d) * (b + d);
    return (n2 * (a + d) - tmp) / (n2 * n2 - tmp);
}

double norm_mutual_info(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes)
{
    double h1 = 0.0, h2 = 0.0, mi = 0.0;
    struct link *link;
    struct graph *g;
    double *p1, *p2;
    double tmp;
    uint32_t i;

    if (!(g = intersection_matrix(labels_true, labels_pred, num_nodes)))
        return NAN;

    if (!(p1 = graph_out_weights(g)))
    {
        free_graph(g);
        return NAN;
    }

    if (!(p2 = graph_in_weights(g)))
    {
        free_graph(g);
        free(p1);
        return NAN;
    }

    for (i = 0; i < g->num_nodes; i++)
    {
        if (p1[i] <= 0.0) continue;
        p1[i] /= num_nodes;
        h1 -= p1[i] * log(p1[i]);
    }

    for (i = 0; i < g->num_nodes; i++)
    {
        if (p2[i] <= 0.0) continue;
        p2[i] /= num_nodes;
        h2 -= p2[i] * log(p2[i]);
    }

    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        if (link->weight <= 0.0) continue;
        tmp = (double)link->weight / num_nodes;
        mi += tmp * log(tmp / (p1[i] * p2[link->index]));
    }

    if (h1 != 0.0 || h2 != 0.0)
        mi = 2.0 * mi / (h1 + h2);
    else
        mi = 1.0;

    free_graph(g);
    free(p1);
    free(p2);
    return mi;
}

void free_labels(uint32_t *labels)
{
    free(labels);
}
