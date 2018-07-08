/*
 * Minimal graph library
 * Label vector manipulation functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#define _GNU_SOURCE
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

int simplify_labels_inplace(uint32_t *labels, uint32_t num_nodes)
{
    uint32_t *lookup;
    uint32_t i, count = 0;

    if (!(lookup = xmalloc(sizeof(*lookup) * num_nodes)))
        return 0;

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
    return 1;
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

void free_labels(uint32_t *labels)
{
    free(labels);
}
