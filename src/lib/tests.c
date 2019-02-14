/*
 * Minimal graph library
 *
 * Copyright (c) 2017-2019 Sebastian Lackner
 */

#include "graph.h"
#include "internal.h"

static void test_xmalloc(void)
{
    int *p;

    p = xmalloc(sizeof(*p));
    assert(p != NULL);
    *p = 13;

    p = xrealloc(p, sizeof(*p) * 2);
    assert(p != NULL);
    assert(*p == 13);
    free(p);
}

static void test_get_edge(void)
{
    struct graph *g = alloc_graph(11, GRAPH_FLAGS_DIRECTED);
    float weight;
    uint32_t i;

    assert(!g->userdata);
    graph_add_edge(g, 0, 2, 1.0);
    graph_add_edge(g, 0, 4, 2.0);
    graph_add_edge(g, 0, 6, 3.0);
    graph_add_edge(g, 0, 8, 4.0);

    for (i = 0; i < g->num_nodes; i++)
    {
        weight = graph_get_edge(g, 0, i);
        if (i < 2 || i > 8 || (i & 1)) assert(weight == 0.0);
        else assert(weight == i / 2.0);
    }

    for (i = g->num_nodes; i < g->num_nodes + 10; i++)
    {
        graph_set_edge(g, 0, i, 1.0);
        graph_set_edge(g, i, 0, 1.0);
        assert(!graph_has_edge(g, 0, i));
        assert(!graph_has_edge(g, i, 0));
        assert(graph_get_edge(g, 0, i) == 0.0);
        assert(graph_get_edge(g, i, 0) == 0.0);
        graph_del_edge(g, 0, i);
        graph_del_edge(g, i, 0);
    }

    free_graph(g);
}

static void test_for_each_edge_directed(void)
{
    struct adjacency empty = {(void *)0xdeadbeef, 0, 0, ~0U};
    struct graph *g = alloc_graph(4, GRAPH_FLAGS_DIRECTED);
    struct link *link;
    double expected;
    uint32_t i;

    graph_add_edge(g, 0, 3, 4.0);
    graph_add_edge(g, 0, 2, 3.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 1, 0, 5.0);
    graph_add_edge(g, 2, 0, 6.0);

    expected = 1.0;
    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        assert(link->weight == expected);
        expected += 1.0;
    }
    assert(expected == 7.0);

    expected = 6.0;
    GRAPH_FOR_EACH_EDGE_REV(g, i, link)
    {
        assert(link->weight == expected);
        expected -= 1.0;
    }
    assert(expected == 0.0);

    GRAPH_FOR_EACH_LINK(&empty, link)
    {
        assert(0);
    }

    expected = 1.0;
    GRAPH_FOR_EACH_LINK(&g->nodes[0], link)
    {
        assert(link->weight == expected);
        expected += 1.0;
    }
    assert(expected == 5.0);

    GRAPH_FOR_EACH_LINK_REV(&empty, link)
    {
        assert(0);
    }

    expected = 4.0;
    GRAPH_FOR_EACH_LINK_REV(&g->nodes[0], link)
    {
        assert(link->weight == expected);
        expected -= 1.0;
    }
    assert(expected == 0.0);

    expected = 0.0;
    GRAPH_FOR_EACH_LINK_ANY(g, &empty, i, link)
    {
        assert(!link);
        assert(i == expected);
        expected += 1.0;
    }
    assert(expected == 4.0);

    expected = 0.0;
    GRAPH_FOR_EACH_LINK_ANY(g, &g->nodes[1], i, link)
    {
        if (i == 0) assert(link->weight == 5.0);
        else assert(!link);
        assert(i == expected);
        expected += 1.0;
    }
    assert(expected == 4.0);

    free_graph(g);
}

static void test_for_each_edge_undirected(void)
{
    struct graph *g = alloc_graph(4, 0);
    struct link *link;
    double expected;
    uint32_t i;

    graph_add_edge(g, 0, 3, 5.0);
    graph_add_edge(g, 0, 2, 3.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 2, 2, 4.0);

    expected = 1.0;
    GRAPH_FOR_EACH_EDGE_UNDIRECTED(g, i, link)
    {
        assert(link->weight == expected);
        expected += 1.0;
    }
    assert(expected == 6.0);

    free_graph(g);
}

static void test_for_each_edge2(void)
{
    struct adjacency empty = {(void *)0xdeadbeef, 0, 0, ~0U};
    struct link *link1, *link2;
    struct graph *g = alloc_graph(4, GRAPH_FLAGS_DIRECTED);
    double expected;
    int i;

    graph_add_edge(g, 0, 1, 1.0);
    graph_add_edge(g, 0, 3, 2.0);
    graph_add_edge(g, 1, 2, 3.0);
    graph_add_edge(g, 1, 3, 4.0);

    GRAPH_FOR_EACH_LINK2(&empty, link1, &empty, link2)
    {
        assert(0);
    }

    expected = 1.0;
    GRAPH_FOR_EACH_LINK2(&g->nodes[0], link1, &empty, link2)
    {
        assert(link1 && link1->weight == expected);
        assert(!link2);
        expected += 1.0;
    }
    assert(expected == 3.0);

    expected = 1.0;
    GRAPH_FOR_EACH_LINK2(&g->nodes[0], link1, &g->nodes[0], link2)
    {
        assert(link1 && link1->weight == expected);
        assert(link1 == link2);
        expected += 1.0;
    }
    assert(expected == 3.0);

    expected = 3.0;
    GRAPH_FOR_EACH_LINK2(&empty, link1, &g->nodes[1], link2)
    {
        assert(link2 && link2->weight == expected);
        assert(!link1);
        expected += 1.0;
    }
    assert(expected == 5.0);

    expected = 3.0;
    GRAPH_FOR_EACH_LINK2(&g->nodes[1], link1, &g->nodes[1], link2)
    {
        assert(link2 && link2->weight == expected);
        assert(link1 == link2);
        expected += 1.0;
    }
    assert(expected == 5.0);

    i = 0;
    GRAPH_FOR_EACH_LINK2(&g->nodes[0], link1, &g->nodes[1], link2)
    {
        switch (i++)
        {
            case 0:
                assert(link1 && link1->index == 1);
                assert(!link2);
                break;
            case 1:
                assert(!link1);
                assert(link2 && link2->index == 2);
                break;
            case 2:
                assert(link1 && link1->index == 3);
                assert(link2 && link2->index == 3);
                break;
            default:
                assert(0);
        }
    }
    assert(i == 3);

    free_graph(g);
}

static void test_for_each_nonexistent_edge(void)
{
    struct adjacency empty = {(void *)0xdeadbeef, 0, 0, ~0U};
    struct graph *g = alloc_graph(4, GRAPH_FLAGS_DIRECTED);
    uint32_t i, expected;

    graph_add_edge(g, 0, 1, 1.0);
    graph_add_edge(g, 0, 3, 2.0);
    graph_add_edge(g, 1, 2, 3.0);
    graph_add_edge(g, 1, 3, 4.0);

    expected = 0;
    GRAPH_FOR_EACH_NONEXISTENT_LINK(g, &empty, i)
    {
        assert(i == expected);
        expected++;
    }
    assert(expected == 4);

    expected = 0;
    GRAPH_FOR_EACH_NONEXISTENT_LINK(g, &g->nodes[0], i)
    {
        assert(i == expected);
        expected += 2;
    }
    assert(expected == 4);

    expected = 0;
    GRAPH_FOR_EACH_NONEXISTENT_LINK(g, &g->nodes[1], i)
    {
        assert(i == expected);
        expected++;
    }
    assert(expected == 2);

    free_graph(g);
}

static void test_assign_graph(void)
{
    struct graph *g = alloc_graph(2, GRAPH_FLAGS_DIRECTED);
    struct graph *h = alloc_graph(2, GRAPH_FLAGS_DIRECTED);

    g->userdata = (void *)0xcafebabe;
    graph_add_edge(g, 0, 0, 10.0);
    graph_add_edge(g, 0, 1, 10.0);
    graph_add_edge(g, 1, 0, 10.0);
    graph_add_edge(g, 1, 1, 10.0);

    h->userdata = (void *)0xdeadbeef;
    graph_add_edge(h, 0, 0, 1.0);
    graph_add_edge(h, 0, 1, 2.0);
    graph_add_edge(h, 1, 0, 3.0);
    graph_add_edge(h, 1, 1, 4.0);

    assign_graph(g, h);

    assert(g->flags == GRAPH_FLAGS_DIRECTED);
    assert(g->userdata == (void *)0xcafebabe);
    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 3.0);
    assert(graph_get_edge(g, 1, 1) == 4.0);

    free_graph(g);
}

static void test_resize_graph(void)
{
    struct graph *g = alloc_graph(2, GRAPH_FLAGS_DIRECTED);
    struct link *link;
    double expected;
    uint32_t i;
    int res;

    g->userdata = (void *)0xdeadbeef;
    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 4.0);
    graph_add_edge(g, 1, 1, 5.0);

    res = resize_graph_inplace(g, 1024);
    assert(res);
    assert(g->userdata == (void *)0xdeadbeef);
    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 4.0);
    assert(graph_get_edge(g, 1, 1) == 5.0);

    graph_add_edge(g, 0, 1023, 3.0);
    graph_add_edge(g, 1, 1023, 6.0);
    graph_add_edge(g, 1023, 0, 7.0);

    expected = 1.0;
    GRAPH_FOR_EACH_EDGE(g, i, link)
    {
        assert(link->weight == expected);
        expected += 1.0;
    }
    assert(expected == 8.0);

    res = resize_graph_inplace(g, 2);
    assert(res);
    assert(g->userdata == (void *)0xdeadbeef);
    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 4.0);
    assert(graph_get_edge(g, 1, 1) == 5.0);
    assert(g->nodes[0].num_links == 2);
    assert(g->nodes[1].num_links == 2);

    res = resize_graph_inplace(g, 0);
    assert(res);
    assert(g->userdata == (void *)0xdeadbeef);
    assert(g->num_nodes == 0);

    free_graph(g);
}

static void test_compress_graph(void)
{
    struct graph *g = alloc_graph(10, GRAPH_FLAGS_DIRECTED);

    graph_add_edge(g, 0, 3, 1.0);
    graph_add_edge(g, 0, 2, 2.0);
    graph_add_edge(g, 0, 1, 3.0);
    assert(g->nodes[0].num_links < g->nodes[0].max_links);

    compress_graph_inplace(g);
    assert(g->nodes[0].num_links == g->nodes[0].max_links);

    free_graph(g);
}

static void test_count_links_directed(void)
{
    struct graph *g = alloc_graph(2, GRAPH_FLAGS_DIRECTED);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);

    assert(graph_count_links(g) == 3);
    assert(graph_count_links_directed(g) == 3);

    free_graph(g);
}

static void test_count_links_undirected(void)
{
    struct graph *g = alloc_graph(2, 0);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);

    assert(graph_count_links(g) == 2);
    assert(graph_count_links_directed(g) == 3);

    free_graph(g);
}

static void test_sum_weights_directed(void)
{
    struct graph *g = alloc_graph(2, GRAPH_FLAGS_DIRECTED);
    uint32_t *degrees;
    double *weights;

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);

    assert(graph_sum_weights(g) == 6.0);

    degrees = graph_out_degrees(g);
    assert(degrees != NULL);
    assert(degrees[0] == 2);
    assert(degrees[1] == 1);
    free(degrees);

    weights = graph_out_weights(g);
    assert(weights != NULL);
    assert(weights[0] == 3.0);
    assert(weights[1] == 3.0);
    free(weights);

    degrees = graph_in_degrees(g);
    assert(degrees != NULL);
    assert(degrees[0] == 2);
    assert(degrees[1] == 1);
    free(degrees);

    weights = graph_in_weights(g);
    assert(weights != NULL);
    assert(weights[0] == 4.0);
    assert(weights[1] == 2.0);
    free(weights);

    free_graph(g);
}

static void test_sum_weights_undirected(void)
{
    struct graph *g = alloc_graph(2, 0);
    uint32_t *degrees;
    double *weights;

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);

    assert(graph_sum_weights(g) == 3.0);

    degrees = graph_out_degrees(g);
    assert(degrees != NULL);
    assert(degrees[0] == 2);
    assert(degrees[1] == 1);
    free(degrees);

    weights = graph_out_weights(g);
    assert(weights != NULL);
    assert(weights[0] == 3.0);
    assert(weights[1] == 2.0);
    free(weights);

    degrees = graph_in_degrees(g);
    assert(degrees != NULL);
    assert(degrees[0] == 2);
    assert(degrees[1] == 1);
    free(degrees);

    weights = graph_in_weights(g);
    assert(weights != NULL);
    assert(weights[0] == 3.0);
    assert(weights[1] == 2.0);
    free(weights);

    free_graph(g);
}

static void test_duplicate_directed(void)
{
    const uint32_t flags = GRAPH_FLAGS_DIRECTED | 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);

    assign_graph(g, duplicate_graph(g));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 3.0);
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
    g = alloc_graph(8, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 0, 2, 3.0);

    assign_graph(g, duplicate_graph(g));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 0, 2) == 3.0);
    graph_add_edge(g, 0, 3, 4.0);
    assert(graph_get_edge(g, 0, 3) == 4.0);

    free_graph(g);
}

static void test_duplicate_undirected(void)
{
    const uint32_t flags = 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);

    assign_graph(g, duplicate_graph(g));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 2.0);
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
    g = alloc_graph(8, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 0, 2, 3.0);

    assign_graph(g, duplicate_graph(g));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 0, 2) == 3.0);
    graph_add_edge(g, 0, 3, 4.0);
    assert(graph_get_edge(g, 0, 3) == 4.0);

    free_graph(g);
}

static void test_transpose_directed(void)
{
    const uint32_t flags = GRAPH_FLAGS_DIRECTED | 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);

    assign_graph(g, transpose_graph(g));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 3.0);
    assert(graph_get_edge(g, 1, 0) == 2.0);
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
}

static void test_transpose_undirected(void)
{
    const uint32_t flags = 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);

    assign_graph(g, transpose_graph(g));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 2.0);
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
}

static void test_invert_directed(void)
{
    const uint32_t flags = GRAPH_FLAGS_DIRECTED | 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);

    assign_graph(g, invert_graph(g, 3.0, 1));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 2.0);
    assert(graph_get_edge(g, 0, 1) == 1.0);
    assert(!graph_has_edge(g, 1, 0));
    assert(graph_get_edge(g, 1, 1) == 3.0);

    assign_graph(g, invert_graph(g, 2.0, 1));
    assert(g->flags == flags);

    assert(!graph_has_edge(g, 0, 0));
    assert(graph_get_edge(g, 0, 1) == 1.0);
    assert(graph_get_edge(g, 1, 0) == 2.0);
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
}

static void test_invert_undirected(void)
{
    const uint32_t flags = 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);

    assign_graph(g, invert_graph(g, 2.0, 1));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(!graph_has_edge(g, 0, 1));
    assert(!graph_has_edge(g, 1, 0));
    assert(graph_get_edge(g, 1, 1) == 2.0);

    assign_graph(g, invert_graph(g, 1.0, 1));
    assert(g->flags == flags);

    assert(!graph_has_edge(g, 0, 0));
    assert(graph_get_edge(g, 0, 1) == 1.0);
    assert(graph_get_edge(g, 1, 0) == 1.0);
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
}

static void test_multiply_const_directed(void)
{
    const uint32_t flags = GRAPH_FLAGS_DIRECTED | 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);

    assign_graph(g, multiply_graph_const(g, 2.0));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 2.0);
    assert(graph_get_edge(g, 0, 1) == 4.0);
    assert(graph_get_edge(g, 1, 0) == 6.0);
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
}

static void test_multiply_const_undirected(void)
{
    const uint32_t flags = 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);

    assign_graph(g, multiply_graph_const(g, 2.0));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 2.0);
    assert(graph_get_edge(g, 0, 1) == 4.0);
    assert(graph_get_edge(g, 1, 0) == 4.0);
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
}

static void test_multiply_elementwise_directed(void)
{
    const uint32_t flags = GRAPH_FLAGS_DIRECTED | 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);

    assign_graph(g, multiply_graph_elementwise(g, g));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 4.0);
    assert(graph_get_edge(g, 1, 0) == 9.0);
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
}

static void test_multiply_elementwise_undirected(void)
{
    const uint32_t flags = 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);

    assign_graph(g, multiply_graph_elementwise(g, g));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 4.0);
    assert(graph_get_edge(g, 1, 0) == 4.0);
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
}

static void test_square_directed(void)
{
    const uint32_t flags = GRAPH_FLAGS_DIRECTED | 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);

    assign_graph(g, square_graph(g));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 7.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 3.0);
    assert(graph_get_edge(g, 1, 1) == 6.0);

    free_graph(g);
}

static void test_square_undirected(void)
{
    const uint32_t flags = 0x00008000;
    struct graph *g = alloc_graph(2, flags);

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);

    assign_graph(g, square_graph(g));
    assert(g->flags == flags);

    assert(graph_get_edge(g, 0, 0) == 5.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 2.0);
    assert(graph_get_edge(g, 1, 1) == 4.0);

    free_graph(g);
}

static void test_filter_directed(void)
{
    struct graph *g = alloc_graph(2, GRAPH_FLAGS_DIRECTED);
    uint32_t labels[2];

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);
    graph_add_edge(g, 1, 1, 4.0);

    labels[0] = 0;
    labels[1] = 1;
    assign_graph(g, filter_graph_labels(g, labels));
    assert(!graph_has_edge(g, 0, 0));
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 3.0);
    assert(!graph_has_edge(g, 1, 1));

    assign_graph(g, filter_graph_weights(g, 2.5, 3.5));
    assert(!graph_has_edge(g, 0, 0));
    assert(!graph_has_edge(g, 0, 1));
    assert(graph_get_edge(g, 1, 0) == 3.0);
    assert(!graph_has_edge(g, 1, 1));

    labels[1] = 0;
    assign_graph(g, filter_graph_labels(g, labels));
    assert(!graph_has_edge(g, 0, 0));
    assert(!graph_has_edge(g, 0, 1));
    assert(!graph_has_edge(g, 1, 0));
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
}

static void test_filter_undirected(void)
{
    struct graph *g = alloc_graph(2, 0);
    uint32_t labels[2];

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 1, 3.0);

    labels[0] = 0;
    labels[1] = 1;
    assign_graph(g, filter_graph_labels(g, labels));
    assert(!graph_has_edge(g, 0, 0));
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 2.0);
    assert(!graph_has_edge(g, 1, 1));

    assign_graph(g, filter_graph_weights(g, 1.5, 2.5));
    assert(!graph_has_edge(g, 0, 0));
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 2.0);
    assert(!graph_has_edge(g, 1, 1));

    labels[1] = 0;
    assign_graph(g, filter_graph_labels(g, labels));
    assert(!graph_has_edge(g, 0, 0));
    assert(!graph_has_edge(g, 0, 1));
    assert(!graph_has_edge(g, 1, 0));
    assert(!graph_has_edge(g, 1, 1));

    free_graph(g);
}

static void test_clamp_directed(void)
{
    struct graph *g = alloc_graph(2, GRAPH_FLAGS_DIRECTED);
    struct graph *h;

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);
    graph_add_edge(g, 1, 1, 4.0);

    h = clamp_graph(g, 2.0, 3.0);
    assert(graph_get_edge(h, 0, 0) == 2.0);
    assert(graph_get_edge(h, 0, 1) == 2.0);
    assert(graph_get_edge(h, 1, 0) == 3.0);
    assert(graph_get_edge(h, 1, 1) == 3.0);
    free_graph(h);

    assert(graph_get_edge(g, 0, 0) == 1.0);
    clamp_graph_inplace(g, 2.0, 3.0);
    assert(graph_get_edge(g, 0, 0) == 2.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 3.0);
    assert(graph_get_edge(g, 1, 1) == 3.0);

    free_graph(g);
}

static void test_clamp_undirected(void)
{
    struct graph *g = alloc_graph(2, 0);
    struct graph *h;

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 1, 3.0);

    h = clamp_graph(g, 1.0, 2.0);
    assert(graph_get_edge(h, 0, 0) == 1.0);
    assert(graph_get_edge(h, 0, 1) == 2.0);
    assert(graph_get_edge(h, 1, 0) == 2.0);
    assert(graph_get_edge(h, 1, 1) == 2.0);
    free_graph(h);

    assert(graph_get_edge(g, 1, 1) == 3.0);
    clamp_graph_inplace(g, 1.0, 2.0);
    assert(graph_get_edge(g, 0, 0) == 1.0);
    assert(graph_get_edge(g, 0, 1) == 2.0);
    assert(graph_get_edge(g, 1, 0) == 2.0);
    assert(graph_get_edge(g, 1, 1) == 2.0);

    free_graph(g);
}

static void test_modularity_directed(void)
{
    struct graph *g = alloc_graph(2, GRAPH_FLAGS_DIRECTED);
    uint32_t labels[2];
    uint32_t count;
    double mod;

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);
    graph_add_edge(g, 1, 0, 3.0);

    labels[0] = 0;
    labels[1] = 0;
    mod = modularity(g, labels);
    assert(mod == 0.0);
    count = count_labels(labels, g->num_nodes);
    assert(count == 1);

    labels[1] = 1;
    mod = modularity(g, labels);
    assert(mod == -1.0 / 3.0);
    count = count_labels(labels, g->num_nodes);
    assert(count == 2);

    labels[0] = 1;
    labels[1] = 0;
    mod = modularity(g, labels);
    assert(mod == -1.0 / 3.0);
    count = count_labels(labels, g->num_nodes);
    assert(count == 2);

    free_graph(g);
}

static void test_modularity_undirected(void)
{
    struct graph *g = alloc_graph(2, 0);
    uint32_t labels[2];
    uint32_t count;
    double mod;

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);

    labels[0] = 0;
    labels[1] = 0;
    mod = modularity(g, labels);
    assert(mod == 0.0);
    count = count_labels(labels, g->num_nodes);
    assert(count == 1);

    labels[1] = 1;
    mod = modularity(g, labels);
    assert(mod == -0.32);
    count = count_labels(labels, g->num_nodes);
    assert(count == 2);

    labels[0] = 1;
    labels[1] = 0;
    mod = modularity(g, labels);
    assert(mod == -0.32);
    count = count_labels(labels, g->num_nodes);
    assert(count == 2);

    free_graph(g);
}

static void test_merge_nodes_add_directed(void)
{
    struct graph *g = alloc_graph(4, GRAPH_FLAGS_DIRECTED);
    uint32_t i, j, c = 0;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            graph_add_edge(g, i, j, c++);

    graph_merge_nodes_add(g, 0, 2);

    /* invalid merges */
    graph_merge_nodes_add(g, 0, 0);
    graph_merge_nodes_add(g, 0, 4);
    graph_merge_nodes_add(g, 4, 0);

    for (i = 0; i < 4; i++)
    {
        assert(!graph_has_edge(g, i, 2));
        assert(!graph_has_edge(g, 2, i));
    }

    assert(graph_get_edge(g, 0, 0) == 20.0);
    assert(graph_get_edge(g, 0, 1) == 10.0);
    assert(graph_get_edge(g, 0, 3) == 14.0);
    assert(graph_get_edge(g, 1, 0) == 10.0);
    assert(graph_get_edge(g, 1, 1) ==  5.0);
    assert(graph_get_edge(g, 1, 3) ==  7.0);
    assert(graph_get_edge(g, 3, 0) == 26.0);
    assert(graph_get_edge(g, 3, 1) == 13.0);
    assert(graph_get_edge(g, 3, 3) == 15.0);

    free_graph(g);
}

static void test_merge_nodes_add_undirected(void)
{
    struct graph *g = alloc_graph(4, 0);
    uint32_t i, j, c = 0;

    for (i = 0; i < 4; i++)
        for (j = 0; j <= i; j++)
            graph_add_edge(g, i, j, c++);

    graph_merge_nodes_add(g, 0, 2);

    /* invalid merges */
    graph_merge_nodes_add(g, 0, 0);
    graph_merge_nodes_add(g, 0, 4);
    graph_merge_nodes_add(g, 4, 0);

    for (i = 0; i < 4; i++)
    {
        assert(!graph_has_edge(g, i, 2));
        assert(!graph_has_edge(g, 2, i));
    }

    assert(graph_get_edge(g, 0, 0) == 11.0);
    assert(graph_get_edge(g, 0, 1) ==  5.0);
    assert(graph_get_edge(g, 0, 3) == 14.0);
    assert(graph_get_edge(g, 1, 0) ==  5.0);
    assert(graph_get_edge(g, 1, 1) ==  2.0);
    assert(graph_get_edge(g, 1, 3) ==  7.0);
    assert(graph_get_edge(g, 3, 0) == 14.0);
    assert(graph_get_edge(g, 3, 1) ==  7.0);
    assert(graph_get_edge(g, 3, 3) ==  9.0);

    free_graph(g);
}

static void test_merge_nodes_min_directed(void)
{
    struct graph *g = alloc_graph(4, GRAPH_FLAGS_DIRECTED);
    uint32_t i, j, c = 0;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            graph_add_edge(g, i, j, c++);

    graph_merge_nodes_min(g, 0, 2);

    /* invalid merges */
    graph_merge_nodes_min(g, 0, 0);
    graph_merge_nodes_min(g, 0, 4);
    graph_merge_nodes_min(g, 4, 0);

    for (i = 0; i < 4; i++)
    {
        assert(!graph_has_edge(g, i, 2));
        assert(!graph_has_edge(g, 2, i));
    }

    assert(graph_get_edge(g, 0, 0) ==  0.0);
    assert(graph_get_edge(g, 0, 1) ==  1.0);
    assert(graph_get_edge(g, 0, 3) ==  3.0);
    assert(graph_get_edge(g, 1, 0) ==  4.0);
    assert(graph_get_edge(g, 1, 1) ==  5.0);
    assert(graph_get_edge(g, 1, 3) ==  7.0);
    assert(graph_get_edge(g, 3, 0) == 12.0);
    assert(graph_get_edge(g, 3, 1) == 13.0);
    assert(graph_get_edge(g, 3, 3) == 15.0);

    free_graph(g);
}

static void test_merge_nodes_min_undirected(void)
{
    struct graph *g = alloc_graph(4, 0);
    uint32_t i, j, c = 0;

    for (i = 0; i < 4; i++)
        for (j = 0; j <= i; j++)
            graph_add_edge(g, i, j, c++);

    graph_merge_nodes_min(g, 0, 2);

    /* invalid merges */
    graph_merge_nodes_min(g, 0, 0);
    graph_merge_nodes_min(g, 0, 4);
    graph_merge_nodes_min(g, 4, 0);

    for (i = 0; i < 4; i++)
    {
        assert(!graph_has_edge(g, i, 2));
        assert(!graph_has_edge(g, 2, i));
    }

    assert(graph_get_edge(g, 0, 0) ==  0.0);
    assert(graph_get_edge(g, 0, 1) ==  1.0);
    assert(graph_get_edge(g, 0, 3) ==  6.0);
    assert(graph_get_edge(g, 1, 0) ==  1.0);
    assert(graph_get_edge(g, 1, 1) ==  2.0);
    assert(graph_get_edge(g, 1, 3) ==  7.0);
    assert(graph_get_edge(g, 3, 0) ==  6.0);
    assert(graph_get_edge(g, 3, 1) ==  7.0);
    assert(graph_get_edge(g, 3, 3) ==  9.0);

    free_graph(g);
}

static void test_merge_nodes_max_directed(void)
{
    struct graph *g = alloc_graph(4, GRAPH_FLAGS_DIRECTED);
    uint32_t i, j, c = 0;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            graph_add_edge(g, i, j, c++);

    graph_merge_nodes_max(g, 0, 2);

    /* invalid merges */
    graph_merge_nodes_max(g, 0, 0);
    graph_merge_nodes_max(g, 0, 4);
    graph_merge_nodes_max(g, 4, 0);

    for (i = 0; i < 4; i++)
    {
        assert(!graph_has_edge(g, i, 2));
        assert(!graph_has_edge(g, 2, i));
    }

    assert(graph_get_edge(g, 0, 0) == 10.0);
    assert(graph_get_edge(g, 0, 1) ==  9.0);
    assert(graph_get_edge(g, 0, 3) == 11.0);
    assert(graph_get_edge(g, 1, 0) ==  6.0);
    assert(graph_get_edge(g, 1, 1) ==  5.0);
    assert(graph_get_edge(g, 1, 3) ==  7.0);
    assert(graph_get_edge(g, 3, 0) == 14.0);
    assert(graph_get_edge(g, 3, 1) == 13.0);
    assert(graph_get_edge(g, 3, 3) == 15.0);

    free_graph(g);
}

static void test_merge_nodes_max_undirected(void)
{
    struct graph *g = alloc_graph(4, 0);
    uint32_t i, j, c = 0;

    for (i = 0; i < 4; i++)
        for (j = 0; j <= i; j++)
            graph_add_edge(g, i, j, c++);

    graph_merge_nodes_max(g, 0, 2);

    /* invalid merges */
    graph_merge_nodes_max(g, 0, 0);
    graph_merge_nodes_max(g, 0, 4);
    graph_merge_nodes_max(g, 4, 0);

    for (i = 0; i < 4; i++)
    {
        assert(!graph_has_edge(g, i, 2));
        assert(!graph_has_edge(g, 2, i));
    }

    assert(graph_get_edge(g, 0, 0) ==  5.0);
    assert(graph_get_edge(g, 0, 1) ==  4.0);
    assert(graph_get_edge(g, 0, 3) ==  8.0);
    assert(graph_get_edge(g, 1, 0) ==  4.0);
    assert(graph_get_edge(g, 1, 1) ==  2.0);
    assert(graph_get_edge(g, 1, 3) ==  7.0);
    assert(graph_get_edge(g, 3, 0) ==  8.0);
    assert(graph_get_edge(g, 3, 1) ==  7.0);
    assert(graph_get_edge(g, 3, 3) ==  9.0);

    free_graph(g);
}

static void test_merge_nodes_custom1(void)  /* single linkage */
{
    struct graph *g = alloc_graph(4, 0);
    uint32_t i, j, c = 0;

    for (i = 0; i < 4; i++)
        for (j = 0; j <= i; j++)
            graph_add_edge(g, i, j, c++);

    graph_merge_nodes_custom(g, 0, 2, 0.5, 0.5, 0.0, -0.5);

    /* invalid merges */
    graph_merge_nodes_custom(g, 0, 0, 0.5, 0.5, 0.0, -0.5);
    graph_merge_nodes_custom(g, 0, 4, 0.5, 0.5, 0.0, -0.5);
    graph_merge_nodes_custom(g, 4, 0, 0.5, 0.5, 0.0, -0.5);

    for (i = 0; i < 4; i++)
    {
        assert(!graph_has_edge(g, i, 2));
        assert(!graph_has_edge(g, 2, i));
    }

    assert(!graph_has_edge(g, 0, 0));
    assert(graph_get_edge(g, 0, 1) ==  1.0);
    assert(graph_get_edge(g, 0, 3) ==  6.0);
    assert(graph_get_edge(g, 1, 0) ==  1.0);
    assert(graph_get_edge(g, 1, 1) ==  2.0);
    assert(graph_get_edge(g, 1, 3) ==  7.0);
    assert(graph_get_edge(g, 3, 0) ==  6.0);
    assert(graph_get_edge(g, 3, 1) ==  7.0);
    assert(graph_get_edge(g, 3, 3) ==  9.0);

    free_graph(g);
}

static void test_merge_nodes_custom2(void)  /* complete linkage */
{
    struct graph *g = alloc_graph(4, 0);
    uint32_t i, j, c = 0;

    for (i = 0; i < 4; i++)
        for (j = 0; j <= i; j++)
            graph_add_edge(g, i, j, c++);

    graph_merge_nodes_custom(g, 0, 2, 0.5, 0.5, 0.0, 0.5);

    /* invalid merges */
    graph_merge_nodes_custom(g, 0, 0, 0.5, 0.5, 0.0, 0.5);
    graph_merge_nodes_custom(g, 0, 4, 0.5, 0.5, 0.0, 0.5);
    graph_merge_nodes_custom(g, 4, 0, 0.5, 0.5, 0.0, 0.5);

    for (i = 0; i < 4; i++)
    {
        assert(!graph_has_edge(g, i, 2));
        assert(!graph_has_edge(g, 2, i));
    }

    assert(!graph_has_edge(g, 0, 0));
    assert(graph_get_edge(g, 0, 1) ==  4.0);
    assert(graph_get_edge(g, 0, 3) ==  8.0);
    assert(graph_get_edge(g, 1, 0) ==  4.0);
    assert(graph_get_edge(g, 1, 1) ==  2.0);
    assert(graph_get_edge(g, 1, 3) ==  7.0);
    assert(graph_get_edge(g, 3, 0) ==  8.0);
    assert(graph_get_edge(g, 3, 1) ==  7.0);
    assert(graph_get_edge(g, 3, 3) ==  9.0);

    free_graph(g);
}

static void test_merge_nodes_custom3(void)  /* median */
{
    struct graph *g = alloc_graph(4, 0);
    uint32_t i, j, c = 0;

    for (i = 0; i < 4; i++)
        for (j = 0; j <= i; j++)
            graph_add_edge(g, i, j, c++);

    graph_merge_nodes_custom(g, 0, 2, 0.5, 0.5, -0.25, 0.0);

    /* invalid merges */
    graph_merge_nodes_custom(g, 0, 0, 0.5, 0.5, -0.25, 0.0);
    graph_merge_nodes_custom(g, 0, 4, 0.5, 0.5, -0.25, 0.0);
    graph_merge_nodes_custom(g, 4, 0, 0.5, 0.5, -0.25, 0.0);

    for (i = 0; i < 4; i++)
    {
        assert(!graph_has_edge(g, i, 2));
        assert(!graph_has_edge(g, 2, i));
    }

    assert(!graph_has_edge(g, 0, 0));
    assert(graph_get_edge(g, 0, 1) ==  1.75);
    assert(graph_get_edge(g, 0, 3) ==  6.25);
    assert(graph_get_edge(g, 1, 0) ==  1.75);
    assert(graph_get_edge(g, 1, 1) ==  2.0);
    assert(graph_get_edge(g, 1, 3) ==  7.0);
    assert(graph_get_edge(g, 3, 0) ==  6.25);
    assert(graph_get_edge(g, 3, 1) ==  7.0);
    assert(graph_get_edge(g, 3, 3) ==  9.0);

    free_graph(g);
}

static void test_antimodularity(void)
{
    struct graph *g = alloc_graph(2, 0);
    uint32_t labels[2];
    double mod;

    graph_add_edge(g, 0, 0, 1.0);
    graph_add_edge(g, 0, 1, 2.0);

    labels[0] = 0;
    labels[1] = 0;
    mod = antimodularity(g, labels);
    assert(mod == 0.25);

    free_graph(g);
}

/*
static void test_sort_adjacency(void)
{
    struct adjacency *adj, *adj2;

    adj2 = sort_adjacency(NULL);
    assert(adj2 == NULL);

    adj = xmalloc(offsetof(struct adjacency, links[0]));
    adj->num_links = 0;
    adj->max_links = 0;
    adj->sorted = 0;

    adj2 = sort_adjacency(adj);
    assert(adj2 == adj);
    assert(adj->sorted);

    free(adj);
    adj = xmalloc(offsetof(struct adjacency, links[2]));
    adj->num_links = 2;
    adj->max_links = 2;
    adj->sorted = 0;
    adj->links[0].index = 1;
    adj->links[0].weight = 1.0;
    adj->links[1].index = 0;
    adj->links[1].weight = 2.0;

    adj2 = sort_adjacency(adj);
    assert(adj2 == adj);
    assert(adj->sorted);
    assert(adj->links[0].index == 0);
    assert(adj->links[0].weight == 2.0);
    assert(adj->links[1].index == 1);
    assert(adj->links[1].weight == 1.0);

    free(adj);
}
*/

static void test_decode_labels(void)
{
    uint32_t indices[3];
    uint32_t *labels;

    indices[0] = 0;
    indices[1] = 0;
    indices[2] = 0;
    labels = decode_labels(indices, 3);
    assert(labels != NULL);
    assert(labels[0] == 0);
    assert(labels[1] == 0);
    assert(labels[2] == 0);
    free(labels);

    indices[0] = 1;
    indices[1] = 0;
    indices[2] = 2;
    labels = decode_labels(indices, 3);
    assert(labels != NULL);
    assert(labels[0] == 0);
    assert(labels[1] == 0);
    assert(labels[2] == 2);
    free(labels);

    indices[0] = 2;
    indices[1] = 3;
    indices[2] = 2;
    labels = decode_labels(indices, 3);
    assert(labels != NULL);
    assert(labels[0] == 0);
    assert(labels[1] == 1);
    assert(labels[2] == 0);
    free(labels);

    indices[0] = 10;
    indices[1] = 2;
    indices[2] = 1;
    labels = decode_labels(indices, 3);
    assert(labels != NULL);
    assert(labels[0] == 0);
    assert(labels[1] == 1);
    assert(labels[2] == 1);
    free(labels);
}

static void test_bfs(void)
{
    struct graph *g = alloc_graph(5, GRAPH_FLAGS_DIRECTED);
    uint32_t len, *lens;
    double dst, *dsts;

    graph_add_edge(g, 0, 1, 1.0);
    graph_add_edge(g, 1, 2, 1.0);
    graph_add_edge(g, 2, 3, 1.0);
    graph_add_edge(g, 3, 4, 1.5);
    graph_add_edge(g, 2, 4, 1.5);

    len = graph_get_distance_count(g, 5, 0);
    assert(len == ~0U);
    dst = graph_get_distance_weight(g, 5, 0);
    assert(dst == INFINITY);

    len = graph_get_distance_count(g, 4, 0);
    assert(len == ~0U);
    dst = graph_get_distance_weight(g, 4, 0);
    assert(dst == INFINITY);

    len = graph_get_distance_count(g, 0, 0);
    assert(len == 0);
    dst = graph_get_distance_weight(g, 0, 0);
    assert(dst == 0.0);

    len = graph_get_distance_count(g, 0, 4);
    assert(len == 3);
    dst = graph_get_distance_weight(g, 0, 4);
    assert(dst == 3.5);

    lens = graph_get_all_distances_count(g, 0, ~0U);
    assert(lens != NULL);
    assert(lens[0] == 0);
    assert(lens[1] == 1);
    assert(lens[2] == 2);
    assert(lens[3] == 3);
    assert(lens[4] == 3);
    free(lens);

    lens = graph_get_all_distances_count(g, 0, 2);
    assert(lens != NULL);
    assert(lens[0] == 0);
    assert(lens[1] == 1);
    assert(lens[2] == 2);
    assert(lens[3] == ~0U);
    assert(lens[4] == ~0U);
    free(lens);

    dsts = graph_get_all_distances_weight(g, 0, INFINITY);
    assert(dsts != NULL);
    assert(dsts[0] == 0.0);
    assert(dsts[1] == 1.0);
    assert(dsts[2] == 2.0);
    assert(dsts[3] == 3.0);
    assert(dsts[4] == 3.5);
    free(dsts);

    dsts = graph_get_all_distances_weight(g, 0, 2.0);
    assert(dsts != NULL);
    assert(dsts[0] == 0.0);
    assert(dsts[1] == 1.0);
    assert(dsts[2] == 2.0);
    assert(dsts[3] == INFINITY);
    assert(dsts[4] == INFINITY);
    free(dsts);

    graph_del_edge(g, 2, 4);

    len = graph_get_distance_count(g, 0, 4);
    assert(len == 4);
    dst = graph_get_distance_weight(g, 0, 4);
    assert(dst == 4.5);

    free_graph(g);
}

static void test_get_connected_components(void)
{
    struct graph *g = alloc_graph(6, 0);
    uint32_t *comp;

    graph_add_edge(g, 0, 1, 1.0);
    graph_add_edge(g, 1, 3, 1.0);
    graph_add_edge(g, 2, 4, 1.0);
    graph_add_edge(g, 4, 5, 1.0);

    comp = graph_get_connected_components(g);
    assert(comp != NULL);
    assert(comp[0] == 0);
    assert(comp[1] == 0);
    assert(comp[2] == 1);
    assert(comp[3] == 0);
    assert(comp[4] == 1);
    assert(comp[5] == 1);
    free(comp);

    free_graph(g);
}

static void test_split_labels(void)
{
    uint32_t *labels1 = xmalloc(sizeof(*labels1) * 8);
    uint32_t *labels2 = xmalloc(sizeof(*labels2) * 8);
    int i;

    for (i = 0; i < 8; i++)
    {
        labels1[i] = i & 4;
        labels2[i] = i & 2;
    }

    split_labels(labels1, labels2, 8);

    for (i = 0; i < 8; i++)
        assert(labels1[i] == (i >> 1));

    free(labels1);
    free(labels2);
}

static int _sort_double(const void *a, const void *b, void *userdata)
{
    const double *da = a, *db = b;
    assert(userdata == (void *)0xdeadbeef);
    return COMPARE(*da, *db);
}

static void test_minheap(void)
{
    static const struct
    {
        int push;
        double value;
    }
    tests[] =
    {
        { 1, 5.0 },
        { 1, 3.0 },
        { 1, 2.0 },
        { 1, 4.0 },
        { 1, 1.0 },
        { 0, 1.0 },
        { 1, 6.0 },
        { 0, 2.0 },
        { 0, 3.0 },
        { 1, 7.0 },
        { 0, 4.0 },
        { 1, 8.0 },
        { 0, 5.0 },
        { 0, 6.0 },
        { 0, 7.0 },
        { 1, 9.0 },
        { 0, 8.0 },
        { 0, 9.0 },
    };
    struct minheap *h;
    double val;
    int i, res;

    h = alloc_minheap(sizeof(double), _sort_double, (void *)0xdeadbeef);
    assert(h != NULL);

    for (i = 0; i < sizeof(tests)/sizeof(tests[0]); i++)
    {
        if (!tests[i].push) res = minheap_pop(h, &val);
        else res = minheap_push(h, &tests[i].value);
        assert(res && (tests[i].push || val == tests[i].value));
    }

    res = minheap_pop(h, &val);
    assert(!res);

    free_minheap(h);
}

int main(void)
{
    test_xmalloc();
    test_get_edge();
    test_for_each_edge_directed();
    test_for_each_edge_undirected();
    test_for_each_edge2();
    test_for_each_nonexistent_edge();
    test_assign_graph();
    test_resize_graph();
    test_compress_graph();
    test_count_links_directed();
    test_count_links_undirected();
    test_sum_weights_directed();
    test_sum_weights_undirected();
    test_duplicate_directed();
    test_duplicate_undirected();
    test_transpose_directed();
    test_transpose_undirected();
    test_invert_directed();
    test_invert_undirected();
    test_multiply_const_directed();
    test_multiply_const_undirected();
    test_multiply_elementwise_directed();
    test_multiply_elementwise_undirected();
    test_square_directed();
    test_square_undirected();
    test_filter_directed();
    test_filter_undirected();
    test_clamp_directed();
    test_clamp_undirected();
    test_modularity_directed();
    test_modularity_undirected();
    test_merge_nodes_add_directed();
    test_merge_nodes_add_undirected();
    test_merge_nodes_min_directed();
    test_merge_nodes_min_undirected();
    test_merge_nodes_max_directed();
    test_merge_nodes_max_undirected();
    test_merge_nodes_custom1();
    test_merge_nodes_custom2();
    test_merge_nodes_custom3();
    test_antimodularity();
    /* test_sort_adjacency(); */
    test_decode_labels();
    test_bfs();
    test_get_connected_components();
    test_split_labels();
    test_minheap();

    fprintf(stderr, "No test failures found\n");
}
