/*
 * Minimal graph library
 *
 * Copyright (c) 2017-2019 Sebastian Lackner
 */

#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <assert.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

struct link
{
    uint32_t index;
    float    weight;
};

struct adjacency
{
    struct link *links;
    uint32_t    num_links;
    uint32_t    max_links;
    uint32_t    hint;
};

#define GRAPH_FLAGS_DIRECTED 0x00000001 /* graph is directed */

struct graph
{
    uint32_t flags;
    uint32_t num_nodes;
    void    *userdata;
    struct adjacency *nodes;
    struct adjacency buffer[0];
};

struct bfs_entry
{
    double   weight;
    uint32_t count;
    uint32_t from;
    uint32_t to;
};

struct _graph_iter
{
    struct link *link;
    struct link *end_link;
    uint32_t    num_nodes;
};

struct _graph_iter2
{
    struct link *link1;
    struct link *link2;
    struct link *end_link1;
    struct link *end_link2;
};

#define _UNIQUE_VARIABLE3(a, b) (a ## b)
#define _UNIQUE_VARIABLE2(a, b) _UNIQUE_VARIABLE3(a, b)
#define _UNIQUE_VARIABLE(a) _UNIQUE_VARIABLE2(a, __COUNTER__)

static inline struct link *__graph_for_each_link(struct adjacency *adj, struct link **link)
{
    *link = &adj->links[0];
    return &adj->links[adj->num_links];
}

static inline struct link *__graph_for_each_link_rev(struct adjacency *adj, struct link **link)
{
    *link = &adj->links[adj->num_links];
    return &adj->links[0];
}

static inline struct _graph_iter2 __graph_for_each_link2(struct adjacency *adj1, struct adjacency *adj2)
{
    struct _graph_iter2 iter;

    iter.end_link1 = __graph_for_each_link(adj1, &iter.link1);
    iter.end_link2 = __graph_for_each_link(adj2, &iter.link2);

    return iter;
}

static inline int __graph_next_link2(struct _graph_iter2 *iter, struct link **link1, struct link **link2)
{
    if (iter->link1 != iter->end_link1)
    {
        if (iter->link2 != iter->end_link2)
        {
            if (iter->link1->index < iter->link2->index)
            {
                *link1 = iter->link1++;
                *link2 = NULL;
            }
            else if (iter->link1->index > iter->link2->index)
            {
                *link1 = NULL;
                *link2 = iter->link2++;
            }
            else
            {
                *link1 = iter->link1++;
                *link2 = iter->link2++;
            }
            return 1;
        }
        else
        {
            *link1 = iter->link1++;
            *link2 = NULL;
            return 1;
        }
    }
    else
    {
        if (iter->link2 != iter->end_link2)
        {
            *link1 = NULL;
            *link2 = iter->link2++;
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

static inline struct _graph_iter __graph_for_each_link_any(const struct graph *g, struct adjacency *adj, uint32_t *index)
{
    struct _graph_iter iter;

    iter.end_link = __graph_for_each_link(adj, &iter.link);
    iter.num_nodes = g->num_nodes;
    *index = 0;

    return iter;
}

static inline int __graph_next_nonexistent_link(struct _graph_iter *iter, uint32_t *index)
{
    for (; *index < iter->num_nodes; (*index)++)
    {
        if (iter->link == iter->end_link ||
            iter->link->index != *index) return 1;
        iter->link++;
    }

    assert(iter->link == iter->end_link);
    return 0;
}

static inline int __graph_next_link_any(struct _graph_iter *iter, uint32_t *index, struct link **link)
{
    if (*index >= iter->num_nodes)
    {
        assert(iter->link == iter->end_link);
        return 0;
    }

    if (iter->link != iter->end_link &&
        iter->link->index == *index)
    {
        *link = iter->link++;
        return 1;
    }

    *link = NULL;
    return 1;
}

#define _GRAPH_FOR_EACH_LINK(_adj, _link, _end_link) \
    for (struct link *(_end_link) = __graph_for_each_link((_adj), &(_link)); \
         (_link) != (_end_link); (_link)++)

#define GRAPH_FOR_EACH_LINK(_adj, _link) \
    _GRAPH_FOR_EACH_LINK((_adj), (_link), _UNIQUE_VARIABLE(__end_link_))

#define _GRAPH_FOR_EACH_LINK_REV(_adj, _link, _end_link) \
    for (struct link *(_end_link) = __graph_for_each_link_rev((_adj), &(_link)); \
         (_link)-- != (_end_link);)

#define GRAPH_FOR_EACH_LINK_REV(_adj, _link) \
    _GRAPH_FOR_EACH_LINK_REV((_adj), (_link), _UNIQUE_VARIABLE(__end_link_))

#define GRAPH_FOR_EACH_EDGE(_g, _i, _link) \
    for ((_i) = 0; (_i) < (_g)->num_nodes; (_i)++) \
        GRAPH_FOR_EACH_LINK(&(_g)->nodes[(_i)], (_link))

#define GRAPH_FOR_EACH_EDGE_REV(_g, _i, _link) \
    for ((_i) = (_g)->num_nodes; (_i)--;) \
        GRAPH_FOR_EACH_LINK_REV(&(_g)->nodes[(_i)], (_link))

#define _GRAPH_FOR_EACH_EDGE_UNDIRECTED(_g, _i, _link, _end_link) \
    for ((_i) = 0; (_i) < (_g)->num_nodes; (_i)++) \
        for (struct link *(_end_link) = __graph_for_each_link(&(_g)->nodes[(_i)], &(_link)); \
             (_link) != (_end_link) && (_link)->index <= (_i); (_link)++)

#define GRAPH_FOR_EACH_EDGE_UNDIRECTED(_g, _i, _link) \
    _GRAPH_FOR_EACH_EDGE_UNDIRECTED((_g), (_i), (_link), _UNIQUE_VARIABLE(__end_link_))

#define _GRAPH_FOR_EACH_LINK2(_adj1, _link1, _adj2, _link2, _iter) \
    for (struct _graph_iter2 (_iter) = __graph_for_each_link2((_adj1), (_adj2)); \
         __graph_next_link2(&(_iter), &(_link1), &(_link2));)

#define GRAPH_FOR_EACH_LINK2(_adj1, _link1, _adj2, _link2) \
    _GRAPH_FOR_EACH_LINK2((_adj1), (_link1), (_adj2), (_link2), _UNIQUE_VARIABLE(__iter_))

#define _GRAPH_FOR_EACH_NONEXISTENT_LINK(_g, _adj, _i, _iter) \
    for (struct _graph_iter (_iter) = __graph_for_each_link_any((_g), (_adj), &(_i)); \
         __graph_next_nonexistent_link(&(_iter), &(_i)); (_i)++)

#define GRAPH_FOR_EACH_NONEXISTENT_LINK(_g, _adj, _i) \
    _GRAPH_FOR_EACH_NONEXISTENT_LINK((_g), (_adj), (_i), _UNIQUE_VARIABLE(__iter_))

#define _GRAPH_FOR_EACH_LINK_ANY(_g, _adj, _i, _link, _iter) \
    for (struct _graph_iter (_iter) = __graph_for_each_link_any((_g), (_adj), &(_i)); \
         __graph_next_link_any(&(_iter), &(_i), &(_link)); (_i)++)

#define GRAPH_FOR_EACH_LINK_ANY(_g, _adj, _i, _link) \
    _GRAPH_FOR_EACH_LINK_ANY((_g), (_adj), (_i), (_link), _UNIQUE_VARIABLE(__iter_))

/* allocate memory and exit on failure */
void *xmalloc(size_t size);
/* reallocate memory and exit on failure */
void *xrealloc(void *ptr, size_t new_size);

/* allocate a new graph structure without any edges */
struct graph *alloc_graph(uint32_t num_nodes, uint32_t flags);
/* load a graph from a file */
struct graph *load_graph(const char *filename, uint32_t flags);
/* duplicate a graph */
struct graph *duplicate_graph(const struct graph *g);
/* compute the transposed graph */
struct graph *transpose_graph(const struct graph *g);
/* compute the inverted graph */
struct graph *invert_graph(const struct graph *g, double max_weight, int self_loops);
/* multiply graph with a constant */
struct graph *multiply_graph_const(const struct graph *g, double constant);
/* multiply graph with a vector */
double *multiply_graph_vector(const struct graph *g, double *vector);
/* multiply graphs elementwise */
struct graph *multiply_graph_elementwise(const struct graph *g, const struct graph *h);
/* add graphs elementwise */
struct graph *add_graph_elementwise(const struct graph *g, const struct graph *h);
/* Compute scalar product of two graphs */
double scalar_product_graph(const struct graph *g, const struct graph *h);
/* compute the graph multiplication */
struct graph *multiply_graph(const struct graph *g, const struct graph *h);
/* compute the squared graph */
struct graph *square_graph(const struct graph *g);
/* compute eigenvector for largest eigenvalue using power iteration */
double *graph_power_iteration(const struct graph *g, uint32_t num_iterations, double *eigenvalue_out);
/* filter a graph */
struct graph *filter_graph_labels(const struct graph *g, uint32_t *labels);
struct graph *filter_graph_weights(const struct graph *g, float min, float max);
/* clamp weights of a graph */
struct graph *clamp_graph(const struct graph *g, float min, float max);
int clamp_graph_inplace(struct graph *g, float min, float max);
/* perform inplace resizing of graph */
int resize_graph_inplace(struct graph *g, uint32_t num_nodes);
/* assign graph and release source */
void assign_graph(struct graph *g, struct graph *src);
/* print a graph adjacency list */
void print_graph(FILE *file, const struct graph *g);
/* try to shrink memory usage of a graph */
void compress_graph_inplace(struct graph *g);
/* deallocate all memory used by a graph */
void free_graph(struct graph *g);

/* check if a graph has a given edge */
int graph_has_edge(const struct graph *g, uint32_t start, uint32_t end);
/* get weight associated with an edge */
float graph_get_edge(const struct graph *g, uint32_t start, uint32_t end);
uint64_t graph_get_edges(const struct graph *g, uint32_t *indices, float *weights, uint64_t max_edges);
void graph_get_weights(const struct graph *g, float *weights);
/* set edges to a given weight */
void graph_set_edge(struct graph *g, uint32_t start, uint32_t end, float weight);
void graph_set_edges(struct graph *g, uint32_t *indices, float *weights, uint64_t num_edges);
/* add weight to a given edge */
void graph_add_edge(struct graph *g, uint32_t start, uint32_t end, float weight);
void graph_add_edges(struct graph *g, uint32_t *indices, float *weights, uint64_t num_edges);
/* add weight to minimum of current edge and given value */
void graph_min_edge(struct graph *g, uint32_t start, uint32_t end, float weight);
void graph_min_edges(struct graph *g, uint32_t *indices, float *weights, uint64_t num_edges);
/* set weight to maximum of current edge and given value */
void graph_max_edge(struct graph *g, uint32_t start, uint32_t end, float weight);
void graph_max_edges(struct graph *g, uint32_t *indices, float *weights, uint64_t num_edges);
/* delete edge from a graph */
void graph_del_edge(struct graph *g, uint32_t start, uint32_t end);
void graph_del_edges(struct graph *g, uint32_t *indices, uint64_t num_edges);

/* merge two nodes in a graph */
void graph_merge_nodes_add(struct graph *g, uint32_t index0, uint32_t index1);
/* merge two nodes in a graph */
void graph_merge_nodes_min(struct graph *g, uint32_t index0, uint32_t index1);
/* merge two nodes in a graph */
void graph_merge_nodes_max(struct graph *g, uint32_t index0, uint32_t index1);
/* merge two nodes in a graph */
void graph_merge_nodes_custom(struct graph *g, uint32_t index0, uint32_t index1,
                              double alpha1, double alpha2, double beta, double gamma);

/* count the number of edges in a graph - COMPLEXITY: O(n) */
uint64_t graph_count_links(const struct graph *g);
/* count the number of edges in a graph (interpreted as directed) - COMPLEXITY: O(n) */
uint64_t graph_count_links_directed(const struct graph *g);
/* count the number of edges on the diagonal */
uint64_t graph_count_links_diagonal(const struct graph *g);
/* sum the weights of all edges in a graph - COMPLEXITY: O(m) */
double graph_sum_weights(const struct graph *g);
/* compute out-degrees for vertices of a graph - COMPLEXITY: O(n) */
uint32_t graph_out_degree(const struct graph *g, uint32_t index);
uint32_t *graph_out_degrees(const struct graph *g);
uint32_t graph_out_edges(const struct graph *g, uint32_t index, uint32_t *indices, float *weights, uint32_t max_edges);
/* compute in-degrees for vertices of a graph - COMPLEXITY: O(m) */
uint32_t graph_in_degree(const struct graph *g, uint32_t index);
uint32_t *graph_in_degrees(const struct graph *g);
/* compute out-weight sums for vertices of a graph - COMPLEXITY: O(m) */
double *graph_out_weights(const struct graph *g);
/* compute in-weight sums for vertices of a graph - COMPLEXITY: O(m) */
double *graph_in_weights(const struct graph *g);
/* compute degree anomalies of a graph */
double *graph_degree_anomalies(const struct graph *g);
/* compute weight anomalies of a graph */
double *graph_weight_anomalies(const struct graph *g);

/* perform BFS on a network */
int graph_bfs(const struct graph *g, uint32_t start, int use_weights,
              int (*callback)(const struct graph *g, struct bfs_entry *entry, void *userdata),
              void *userdata);

/* get number of hops on the path from 'start' to 'end' */
uint32_t graph_get_distance_count(const struct graph *g, uint32_t start, uint32_t end);
/* get sum of weights on the path from 'start' to 'end' */
double graph_get_distance_weight(const struct graph *g, uint32_t start, uint32_t end);
/* get number of hops to all other nodes */
uint32_t *graph_get_all_distances_count(const struct graph *g, uint32_t start, uint32_t max_count);
/* get sum of weights to all other nodes */
double *graph_get_all_distances_weight(const struct graph *g, uint32_t start, double max_weight);
/* return a graph with distance information of the original graph */
struct graph *graph_get_all_distances_graph(const struct graph *g, int use_weights);
/* get the connected components of a graph */
uint32_t *graph_get_connected_components(const struct graph *g);

/* allocate a new label structure */
uint32_t *alloc_labels(uint32_t num_nodes);
/* decode labels from indicies */
uint32_t *decode_labels(const uint32_t *indices, uint32_t num_nodes);
/* load labels from a file */
uint32_t *load_labels(const char *filename, uint32_t num_nodes);
/* print a list of labels */
void print_labels(FILE *file, uint32_t *labels, uint32_t num_nodes);
/* count number of unique labels */
uint32_t count_labels(const uint32_t *labels, uint32_t num_nodes);
/* simplify labels (replace with small integer numbers) */
uint32_t *simplify_labels(const uint32_t *labels, uint32_t num_nodes, uint32_t *count_out);
uint32_t simplify_labels_inplace(uint32_t *labels, uint32_t num_nodes);
/* split labels (e.g. based on connected components) */
void split_labels(uint32_t *labels1, const uint32_t *labels2, uint32_t num_nodes);
/* compute the intersection matrix using two label vectors */
struct graph *intersection_matrix(const uint32_t *labels1, const uint32_t *labels2, uint32_t num_nodes);
/* compute the confusion matrix using two label vectors */
int confusion_matrix(double *a_out, double *b_out, double *c_out, double *d_out,
                     const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes);
int confusion_matrix_comp(double *a_out, double *b_out, double *c_out, double *d_out,
                          const uint32_t *labels_true, const uint32_t *labels_pred, const struct graph *g);
/* compute pair counting measures */
double precision(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes);
double recall(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes);
double rand_index(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes);
double fowlkes_mallows(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes);
double jaccard(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes);
double f1_measure(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes);
double adj_rand_index(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes);
double adj_rand_index_comp(const uint32_t *labels_true, const uint32_t *labels_pred, const struct graph *g);
/* compute the normalized mutual information */
double norm_mutual_info(const uint32_t *labels_true, const uint32_t *labels_pred, uint32_t num_nodes);
/* deallocate memory for labels */
void free_labels(uint32_t *labels);
/* compute modularity of a graph - COMPLEXITY: O(m + n log(n)) */
double modularity(const struct graph *g, const uint32_t *labels);
/* decode labels and compute modularity */
double modularity_decode(const struct graph *g, const uint32_t *indices);
/* compute anti-modularity of a graph */
double antimodularity(const struct graph *g, const uint32_t *labels);
/* decode labels and compute antimodularity */
double antimodularity_decode(const struct graph *g, const uint32_t *indices);

/* allocate a new minheap structure */
struct minheap *alloc_minheap(size_t size, int (*compar)(const void *, const void *, void *), void *userdata);
/* deallocate memory for minheap structure */
void free_minheap(struct minheap *h);
/* insert element into minheap structure */
int minheap_push(struct minheap *h, const void *element);
/* heapify operation, only for internal usage */
void minheap_heapify(struct minheap *h, size_t i);
/* remove element from the minheap structure */
int minheap_pop(struct minheap *h, void *element);

/* resort a list after changing a single element */
size_t resort_r(void *base, size_t nmemb, size_t size, size_t index,
                int (*compar)(const void *, const void *, void *), void *userdata);

/* set or replace a progress message */
void progress(const char *format, ...) __attribute__((format (printf,1,2)));

/* label propagation algorithm */
uint32_t *laprop(struct graph *g);

/* agglomerative hierarchical clustering */
#define NEWMAN_FLAGS_MINIMIZE       0x00000001 /* maximize result */
#define NEWMAN_FLAGS_FAST           0x00000002
#define NEWMAN_FLAGS_ANTIMODULARITY 0x00000004 /* maintain lookup tables for sorted weights */
#define NEWMAN_FLAGS_VERBOSE        0x00000008 /* verbose output */
#define NEWMAN_FLAGS_STRICT         0x00000010 /* enable assertions */
uint32_t *newman(struct graph *g, uint32_t flags, uint32_t num_clusters);

/* centrality measures */
double *graph_closeness_centrality(const struct graph *g);
double *graph_harmonic_centrality(const struct graph *g, uint32_t max_count);
double *graph_harmonic_centrality_fast(const struct graph *g, uint32_t max_count);
double graph_transitivity(const struct graph *g);
double *graph_clustering_coefficient(const struct graph *g, int use_weights);
double graph_reciprocity(const struct graph *g);

#endif /* _GRAPH_H_ */
