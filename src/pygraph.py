#!/usr/bin/python2
from ctypes import cdll
from ctypes import cast
from ctypes import c_int
from ctypes import c_uint
from ctypes import c_uint64
from ctypes import c_float
from ctypes import c_void_p
from ctypes import c_double
from ctypes import c_char_p
from ctypes import Structure
from ctypes import POINTER
from ctypes import CFUNCTYPE
from ctypes.util import find_library
import numpy as np
import numpy.ctypeslib as npc
import sys
import os

filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib/libgraph.so")
lib = cdll.LoadLibrary(filename)
libc = cdll.LoadLibrary(find_library('c'))

class c_graph(Structure):
    _fields_ = [("flags", c_uint), ("num_nodes", c_uint)]

class c_bfs_entry(Structure):
    _fields_ = [("weight", c_double), ("count", c_uint), ("edge_from", c_uint), ("edge_to", c_uint)]

# Hacky: we need optional ndpointer parameters at some places.
def or_null(t):
    class wrap:
        def from_param(cls, obj):
            if obj is None: return None
            return t.from_param(obj)
    return wrap()

c_graph_p = POINTER(c_graph)
c_bfs_entry_p = POINTER(c_bfs_entry)
c_bfs_callback_p = CFUNCTYPE(c_int, c_graph_p, c_bfs_entry_p, c_void_p)

lib.alloc_graph.argtypes = (c_uint, c_uint)
lib.alloc_graph.restype = c_graph_p

lib.load_graph.argtypes = (c_char_p, c_uint)
lib.load_graph.restype = c_graph_p

lib.duplicate_graph.argtypes = (c_graph_p,)
lib.duplicate_graph.restype = c_graph_p

lib.transpose_graph.argtypes = (c_graph_p,)
lib.transpose_graph.restype = c_graph_p

lib.invert_graph.argtypes = (c_graph_p, c_double, c_int)
lib.invert_graph.restype = c_graph_p

lib.multiply_graph_const.argtypes = (c_graph_p, c_double)
lib.multiply_graph_const.restype = c_graph_p

lib.multiply_graph_elementwise.argtypes = (c_graph_p, c_graph_p)
lib.multiply_graph_elementwise.restype = c_graph_p

lib.scalar_product_graph.argtypes = (c_graph_p, c_graph_p)
lib.scalar_product_graph.restype = c_double

lib.multiply_graph.argtypes = (c_graph_p, c_graph_p)
lib.multiply_graph.restype = c_graph_p

lib.square_graph.argtypes = (c_graph_p,)
lib.square_graph.restype = c_graph_p

lib.filter_graph_labels.argtypes = (c_graph_p, npc.ndpointer(dtype=np.uint32))
lib.filter_graph_labels.restype = c_graph_p

lib.filter_graph_weights.argtypes = (c_graph_p, c_float, c_float)
lib.filter_graph_weights.restype = c_graph_p

lib.clamp_graph.argtypes = (c_graph_p, c_float, c_float)
lib.clamp_graph.restype = c_graph_p

lib.resize_graph_inplace.argtypes = (c_graph_p, c_uint)
lib.resize_graph_inplace.restype = c_graph_p

lib.compress_graph_inline.argtypes = (c_graph_p,)

lib.free_graph.argtypes = (c_graph_p,)

lib.graph_has_edge.argtypes = (c_graph_p, c_uint, c_uint)
lib.graph_has_edge.restype = c_int

lib.graph_get_edge.argtypes = (c_graph_p, c_uint, c_uint)
lib.graph_get_edge.restype = c_float

lib.graph_get_edges.argtypes = (c_graph_p, or_null(npc.ndpointer(dtype=np.uint32)), or_null(npc.ndpointer(dtype=np.float32)), c_uint64)
lib.graph_get_edges.restype = c_uint64

lib.graph_get_weights.argtypes = (c_graph_p, npc.ndpointer(dtype=np.float32))

lib.graph_set_edge.argtypes = (c_graph_p, c_uint, c_uint, c_float)
lib.graph_set_edges.argtypes = (c_graph_p, npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.float32), c_uint64)

lib.graph_add_edge.argtypes = (c_graph_p, c_uint, c_uint, c_float)
lib.graph_add_edges.argtypes = (c_graph_p, npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.float32), c_uint64)

lib.graph_min_edge.argtypes = (c_graph_p, c_uint, c_uint, c_float)
lib.graph_min_edges.argtypes = (c_graph_p, npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.float32), c_uint64)

lib.graph_max_edge.argtypes = (c_graph_p, c_uint, c_uint, c_float)
lib.graph_max_edges.argtypes = (c_graph_p, npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.float32), c_uint64)

lib.graph_del_edge.argtypes = (c_graph_p, c_uint, c_uint)
lib.graph_del_edges.argtypes = (c_graph_p, npc.ndpointer(dtype=np.uint32), c_uint64)

lib.graph_merge_nodes_add.argtypes = (c_graph_p, c_uint, c_uint)

lib.graph_merge_nodes_min.argtypes = (c_graph_p, c_uint, c_uint)

lib.graph_merge_nodes_max.argtypes = (c_graph_p, c_uint, c_uint)

lib.graph_merge_nodes_custom.argtypes = (c_graph_p, c_uint, c_uint, c_double, c_double, c_double, c_double)

lib.graph_count_links.argtypes = (c_graph_p,)
lib.graph_count_links.restype = c_uint64

lib.graph_count_links_directed.argtypes = (c_graph_p,)
lib.graph_count_links_directed.restype = c_uint64

lib.graph_sum_weights.argtypes = (c_graph_p,)
lib.graph_sum_weights.restype = c_double

lib.graph_out_degrees.argtypes = (c_graph_p,)
lib.graph_out_degrees.restype = POINTER(c_uint)

lib.graph_in_degrees.argtypes = (c_graph_p,)
lib.graph_in_degrees.restype = POINTER(c_uint)

lib.graph_out_weights.argtypes = (c_graph_p,)
lib.graph_out_weights.restype = POINTER(c_double)

lib.graph_in_weights.argtypes = (c_graph_p,)
lib.graph_in_weights.restype = POINTER(c_double)

lib.graph_bfs.argtypes = (c_graph_p, c_uint, c_int, c_bfs_callback_p, c_void_p)
lib.graph_bfs.restype = c_int

lib.graph_get_distance_count.argtypes = (c_graph_p, c_uint, c_uint)
lib.graph_get_distance_count.restype = c_uint

lib.graph_get_distance_weight.argtypes = (c_graph_p, c_uint, c_uint)
lib.graph_get_distance_weight.restype = c_double

lib.graph_get_all_distances_count.argtypes = (c_graph_p, c_uint, c_uint)
lib.graph_get_all_distances_count.restype = POINTER(c_uint)

lib.graph_get_all_distances_weight.argtypes = (c_graph_p, c_uint, c_double)
lib.graph_get_all_distances_weight.restype = POINTER(c_double)

lib.graph_get_all_distances_graph.argtypes = (c_graph_p, c_int)
lib.graph_get_all_distances_graph.restype = c_graph_p

lib.graph_get_connected_components.argtypes = (c_graph_p,)
lib.graph_get_connected_components.restype = POINTER(c_uint)

lib.decode_labels.argtypes = (npc.ndpointer(dtype=np.uint32), c_uint)
lib.decode_labels.restype = POINTER(c_uint)

lib.count_labels.argtypes = (npc.ndpointer(dtype=np.uint32), c_uint)
lib.count_labels.restype = c_uint

lib.simplify_labels_inplace.argtypes = (npc.ndpointer(dtype=np.uint32), c_uint)
lib.simplify_labels_inplace.restype = c_uint

lib.split_labels.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)

lib.intersection_matrix.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)
lib.intersection_matrix.restype = c_graph_p

lib.confusion_matrix.argtypes = (POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                 npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)
lib.confusion_matrix.restype = c_int

lib.precision.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)
lib.precision.restype = c_double

lib.recall.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)
lib.recall.restype = c_double

lib.rand_index.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)
lib.rand_index.restype = c_double

lib.fowlkes_mallows.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)
lib.fowlkes_mallows.restype = c_double

lib.jaccard.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)
lib.jaccard.restype = c_double

lib.f1_measure.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)
lib.f1_measure.restype = c_double

lib.adj_rand_index.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)
lib.adj_rand_index.restype = c_double

lib.adj_rand_index_comp.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_graph_p)
lib.adj_rand_index_comp.restype = c_double

lib.norm_mutual_info.argtypes = (npc.ndpointer(dtype=np.uint32), npc.ndpointer(dtype=np.uint32), c_uint)
lib.norm_mutual_info.restype = c_double

lib.free_labels.argtypes = (POINTER(c_uint),)

lib.modularity.argtypes = (c_graph_p, npc.ndpointer(dtype=np.uint32))
lib.modularity.restype = c_double

lib.modularity_decode.argtypes = (c_graph_p, npc.ndpointer(dtype=np.uint32))
lib.modularity_decode.restype = c_double

lib.antimodularity.argtypes = (c_graph_p, npc.ndpointer(dtype=np.uint32))
lib.antimodularity.restype = c_double

lib.antimodularity_decode.argtypes = (c_graph_p, npc.ndpointer(dtype=np.uint32))
lib.antimodularity_decode.restype = c_double

lib.newman.argtypes = (c_graph_p, c_uint, c_uint)
lib.newman.restype = POINTER(c_uint)

lib.graph_closeness_centrality.argtypes = (c_graph_p,)
lib.graph_closeness_centrality.restype = POINTER(c_double)

lib.graph_harmonic_centrality.argtypes = (c_graph_p, c_uint)
lib.graph_harmonic_centrality.restype = POINTER(c_double)

lib.graph_harmonic_centrality_fast.argtypes = (c_graph_p, c_uint)
lib.graph_harmonic_centrality_fast.restype = POINTER(c_double)

lib.graph_transitivity.argtypes = (c_graph_p,)
lib.graph_transitivity.restype = c_double

lib.graph_clustering_coefficient.argtypes = (c_graph_p, c_int)
lib.graph_clustering_coefficient.restype = POINTER(c_double)

lib.graph_reciprocity.argtypes = (c_graph_p,)
lib.graph_reciprocity.restype = c_double

libc.free.argtypes = (c_void_p,)

GRAPH_FLAGS_DIRECTED = 0x00000001
GRAPH_FLAGS_UNSORTED = 0x00000004

NEWMAN_FLAGS_MINIMIZE       = 0x00000001 # maximize result
NEWMAN_FLAGS_FAST           = 0x00000002
NEWMAN_FLAGS_ANTIMODULARITY = 0x00000004 # maintain lookup tables for sorted weights
NEWMAN_FLAGS_VERBOSE        = 0x00000008 # verbose output
NEWMAN_FLAGS_STRICT         = 0x00000010 # enable assertions

class Graph(object):
    def __init__(self, num_nodes=None, directed=False, sorted=True, obj=None):
        if obj is None:
            flags = 0
            flags |= (GRAPH_FLAGS_DIRECTED if directed else 0)
            flags |= (GRAPH_FLAGS_UNSORTED if not sorted else 0)
            obj = lib.alloc_graph(num_nodes, flags)

        self.obj = obj
        if not obj:
            raise ValueError

    def __del__(self):
        if lib is None:
            return
        if self.obj:
            lib.free_graph(self.obj)

    @property
    def num_nodes(self):
        return self.obj.contents.num_nodes

    @property
    def flags(self):
        return self.obj.contents.flags

    @flags.setter
    def flags(self, value):
        self.obj.contents.flags = value

    @property
    def directed(self):
        return (self.obj.contents.flags & GRAPH_FLAGS_DIRECTED) != 0

    @staticmethod
    def load_graph(filename, directed=False, sorted=True):
        flags = 0
        flags |= (GRAPH_FLAGS_DIRECTED if directed else 0)
        flags |= (GRAPH_FLAGS_UNSORTED if not sorted else 0)
        obj = lib.load_graph(filename, flags)
        if not obj:
            raise IOError
        return Graph(obj=obj)

    def copy(self):
        return Graph(obj=lib.duplicate_graph(self.obj))

    def transpose(self):
        return Graph(obj=lib.transpose_graph(self.obj))

    def invert(self, max_weight=1.0, self_loops=True):
        return Graph(obj=lib.invert_graph(self.obj, max_weight, self_loops))

    def multiply_const(self, constant):
        return Graph(obj=lib.multiply_graph_const(self.obj, constant))

    def multiply_elementwise(self, other):
        if self.num_nodes != other.num_nodes:
            raise ValueError("graphs have different number of nodes")

        return Graph(obj=lib.multiply_graph_elementwise(self.obj, other.obj))

    def scalar_product(self, other):
        if self.num_nodes != other.num_nodes:
            raise ValueError("graphs have different number of nodes")

        return lib.scalar_product_graph(self.obj, other.obj)

    def multiply(self, other):
        if self.num_nodes != other.num_nodes:
            raise ValueError("graphs have different number of nodes")

        return Graph(obj=lib.multiply_graph(self.obj, other.obj))

    def square(self):
        return Graph(obj=lib.square_graph(self.obj))

    def filter_labels(self, labels):
        if len(labels.shape) != 1 or labels.shape[0] != self.num_nodes:
            raise ValueError("labels array does not have correct dimensions")

        labels = np.asarray(labels, dtype=np.uint32)
        return Graph(obj=lib.filter_graph_labels(self.obj, labels))

    def filter_weights(self, min, max):
        return Graph(obj=lib.filter_graph_weights(self.obj, min, max))

    def clamp(self, min, max):
        return Graph(obj=lib.clamp_graph(self.obj, min, max))

    def resize(self, num_nodes):
        if not lib.resize_graph_inplace(self.obj, num_nodes):
            raise MemoryError

    def dump(self):
        raise NotImplementedError()

    def compress(self):
        lib.compress_graph_inline(self.obj)

    def __getitem__(self, pos):
        start, end = pos
        return lib.graph_get_edge(self.obj, start, end)

    def __setitem__(self, pos, value):
        start, end = pos
        lib.graph_set_edge(self.obj, start, end, value)

    def set_edges(self, edges, weights):
        if len(edges.shape) != 2 or edges.shape[1] != 2:
            raise ValueError("edges array does not have correct dimensions")
        if len(weights.shape) != 1:
            raise ValueError("weights array does not have correct dimensions")
        if edges.shape[0] != weights.shape[0]:
            raise ValueError("edges/weights arrays have different length")

        edges = np.asarray(edges, dtype=np.uint32)
        weights = np.asarray(weights, dtype=np.float32)
        lib.graph_set_edges(self.obj, edges, weights, edges.shape[0])

    def __delitem__(self, pos):
        start, end = pos
        lib.graph_del_edge(self.obj, start, end)

    def del_edges(self, edges):
        if len(edges.shape) != 2 or edges.shape[1] != 2:
            raise ValueError("edges array does not have correct dimensions")

        edges = np.asarray(edges, dtype=np.uint32)
        lib.graph_del_edges(self.obj, edges, edges.shape[0])

    def has_edge(self, pos):
        start, end = pos
        return lib.graph_has_edge(self.obj, start, end)

    def add_edge(self, pos, value):
        start, end = pos
        lib.graph_add_edge(self.obj, start, end, value)

    def add_edges(self, edges, weights):
        if len(edges.shape) != 2 or edges.shape[1] != 2:
            raise ValueError("edges array does not have correct dimensions")
        if len(weights.shape) != 1:
            raise ValueError("weights array does not have correct dimensions")
        if edges.shape[0] != weights.shape[0]:
            raise ValueError("edges/weights arrays have different length")

        edges = np.asarray(edges, dtype=np.uint32)
        weights = np.asarray(weights, dtype=np.float32)
        lib.graph_add_edges(self.obj, edges, weights, edges.shape[0])

    def min_edges(self, edges, weights):
        if len(edges.shape) != 2 or edges.shape[1] != 2:
            raise ValueError("edges array does not have correct dimensions")
        if len(weights.shape) != 1:
            raise ValueError("weights array does not have correct dimensions")
        if edges.shape[0] != weights.shape[0]:
            raise ValueError("edges/weights arrays have different length")

        edges = np.asarray(edges, dtype=np.uint32)
        weights = np.asarray(weights, dtype=np.float32)
        lib.graph_min_edges(self.obj, edges, weights, edges.shape[0])

    def max_edges(self, edges, weights):
        if len(edges.shape) != 2 or edges.shape[1] != 2:
            raise ValueError("edges array does not have correct dimensions")
        if len(weights.shape) != 1:
            raise ValueError("weights array does not have correct dimensions")
        if edges.shape[0] != weights.shape[0]:
            raise ValueError("edges/weights arrays have different length")

        edges = np.asarray(edges, dtype=np.uint32)
        weights = np.asarray(weights, dtype=np.float32)
        lib.graph_max_edges(self.obj, edges, weights, edges.shape[0])

    # FIXME: Would be better to unify edges() and weights().
    def edges(self):
        max_links = self.count_links()
        edges = np.empty(shape=(max_links, 2), dtype=np.uint32, order='C')
        num_links = lib.graph_get_edges(self.obj, edges, None, max_links)
        assert num_links == max_links
        return edges

    def weights(self):
        max_links = self.count_links()
        weights = np.empty(shape=(max_links,), dtype=np.float32, order='C')
        num_links = lib.graph_get_edges(self.obj, None, weights, max_links)
        assert num_links == max_links
        return weights

    def matrix(self):
        matrix = np.empty(shape=(self.num_nodes, self.num_nodes), dtype=np.float32, order='C')
        lib.graph_get_weights(self.obj, matrix)
        return matrix

    def merge_nodes_add(self, index1, index2):
        lib.graph_merge_nodes_add(self.obj, index1, index2)

    def merge_nodes_min(self, index1, index2):
        lib.graph_merge_nodes_min(self.obj, index1, index2)

    def merge_nodes_max(self, index1, index2):
        lib.graph_merge_nodes_max(self.obj, index1, index2)

    def merge_nodes_custom(self, index1, index2, alpha1=0.0, alpha2=0.0, beta=0.0, gamma=0.0):
        lib.graph_merge_nodes_custom(self.obj, index1, index2, alpha1, alpha2, beta, gamma)

    def count_links(self, directed=None):
        if directed:
            return lib.graph_count_links_directed(self.obj)
        else:
            return lib.graph_count_links(self.obj)

    def sum_weights(self):
        return lib.graph_sum_weights(self.obj)

    def get_out_degrees(self):
        degrees_p = lib.graph_out_degrees(self.obj)
        degrees = npc.as_array(degrees_p, shape=(self.num_nodes,)).copy()
        libc.free(degrees_p)
        return degrees

    def get_in_degrees(self):
        degrees_p = lib.graph_in_degrees(self.obj)
        degrees = npc.as_array(degrees_p, shape=(self.num_nodes,)).copy()
        libc.free(degrees_p)
        return degrees

    def get_out_weights(self):
        degrees_p = lib.graph_out_weights(self.obj)
        degrees = npc.as_array(degrees_p, shape=(self.num_nodes,)).copy()
        libc.free(degrees_p)
        return degrees

    def get_in_weights(self):
        degrees_p = lib.graph_in_weights(self.obj)
        degrees = npc.as_array(degrees_p, shape=(self.num_nodes,)).copy()
        libc.free(degrees_p)
        return degrees

    def bfs_count(self, start, max_count=0xffffffff):
        result = []

        def wrapper(graph, entry, userdata):
            if entry.contents.count > max_count:
                return 1

            entry = entry.contents
            result.append((entry.weight, entry.count, entry.edge_from, entry.edge_to))
            return 0

        lib.graph_bfs(self.obj, start, 0, c_bfs_callback_p(wrapper), None)
        return result

    def bfs_weight(self, start, max_weight=np.inf):
        result = []

        def wrapper(graph, entry, userdata):
            if entry.contents.weight > max_weight:
                return 1

            entry = entry.contents
            result.append((entry.weight, entry.count, entry.edge_from, entry.edge_to))
            return 0

        lib.graph_bfs(self.obj, start, 1, c_bfs_callback_p(wrapper), None)
        return result

    def get_distance_count(self, start, end):
        return lib.graph_get_distance_count(self.obj, start, end)

    def get_distance_weight(self, start, end):
        return lib.graph_get_distance_weight(self.obj, start, end)

    def get_all_distances_count(self, start, max_count=0xffffffff):
        counts_p = lib.graph_get_all_distances_count(self.obj, start, max_count)
        count = npc.as_array(counts_p, shape=(self.num_nodes,)).copy()
        libc.free(counts_p)
        return count

    def get_all_distances_weight(self, start, max_weight=np.inf):
        weights_p = lib.graph_get_all_distances_weight(self.obj, start, max_weight)
        weights = npc.as_array(weights_p, shape=(self.num_nodes,)).copy()
        libc.free(weights_p)
        return weights

    def get_all_distances_graph(self, use_weights=False):
        return Graph(obj=lib.graph_get_all_distances_graph(self.obj, use_weights))

    def get_connected_components(self):
        components_p = lib.graph_get_connected_components(self.obj)
        components = npc.as_array(components_p, shape=(self.num_nodes,)).copy()
        libc.free(components_p)
        return components

    def modularity(self, labels):
        if len(labels.shape) != 1 or labels.shape[0] != self.num_nodes:
            raise ValueError("labels array does not have correct dimensions")

        labels = np.asarray(labels, dtype=np.uint32)
        return lib.modularity(self.obj, labels)

    def modularity_decode(self, indices):
        if len(indices.shape) != 1 or indices.shape[0] != self.num_nodes:
            raise ValueError("indices array does not have correct dimensions")

        indices = np.asarray(indices, dtype=np.uint32)
        return lib.modularity_decode(self.obj, indices)

    def antimodularity(self, labels):
        if len(labels.shape) != 1 or labels.shape[0] != self.num_nodes:
            raise ValueError("labels array does not have correct dimensions")

        labels = np.asarray(labels, dtype=np.uint32)
        return lib.antimodularity(self.obj, labels)

    def antimodularity_decode(self, indices):
        if len(indices.shape) != 1 or indices.shape[0] != self.num_nodes:
            raise ValueError("indices array does not have correct dimensions")

        indices = np.asarray(indices, dtype=np.uint32)
        return lib.antimodularity_decode(self.obj, indices)

    def newman(self, flags, num_clusters):
        temp = self.copy() # Don't modify user graph
        labels_p = lib.newman(temp.obj, flags, num_clusters)
        labels = npc.as_array(labels_p, shape=(self.num_nodes,)).copy()
        lib.free_labels(labels_p)
        return labels

    def closeness_centrality(self):
        centrality_p = lib.graph_closeness_centrality(self.obj)
        centrality = npc.as_array(centrality_p, shape=(self.num_nodes,)).copy()
        libc.free(centrality_p)
        return centrality

    def harmonic_centrality(self, max_count=0xffffffff, fast=False):
        if fast and max_count in [0, 1, 2, 3]:
            centrality_p = lib.graph_harmonic_centrality_fast(self.obj, max_count)
        else:
            centrality_p = lib.graph_harmonic_centrality(self.obj, max_count)
        centrality = npc.as_array(centrality_p, shape=(self.num_nodes,)).copy()
        libc.free(centrality_p)
        return centrality

    def transitivity(self):
        assert not self.directed
        return lib.graph_transitivity(self.obj)

    def clustering_coefficient(self, use_weights=False):
        clustering_p = lib.graph_clustering_coefficient(self.obj, use_weights)
        clustering = npc.as_array(clustering_p, shape=(self.num_nodes,)).copy()
        libc.free(clustering_p)
        return clustering

    def reciprocity(self):
        assert self.directed
        return lib.graph_reciprocity(self.obj)

    def adj_rand_index_comp(self, indices_true, indices_pred):
        if len(indices_true.shape) != 1 or len(indices_pred.shape) != 1:
            raise ValueError("indices array does not have correct dimensions")
        if indices_true.shape[0] != indices_pred.shape[0]:
            raise ValueError("indices arrays have different length")

        indices_true = np.asarray(indices_true, dtype=np.uint32)
        indices_pred = np.asarray(indices_pred, dtype=np.uint32)
        return lib.adj_rand_index_comp(indices_true, indices_pred, self.obj)

def decode_labels(indices):
    if len(indices.shape) != 1:
        raise ValueError("indices array does not have correct dimensions")

    indices = np.asarray(indices, dtype=np.uint32)
    labels_p = lib.decode_labels(indices, indices.shape[0])
    labels = npc.as_array(labels_p, shape=indices.shape).copy()
    lib.free_labels(labels_p)
    return labels

def count_labels(labels):
    if len(labels.shape) != 1:
        raise ValueError("labels array does not have correct dimensions")

    labels = np.asarray(labels, dtype=np.uint32)
    count = lib.count_labels(labels, labels.shape[0])
    if count == 0xffffffff:
        raise MemoryError

    return count

def simplify_labels(indices):
    if len(indices.shape) != 1:
        raise ValueError("indices array does not have correct dimensions")

    indices = np.asarray(indices.copy(), dtype=np.uint32)
    count = lib.simplify_labels_inplace(indices, indices.shape[0])
    if count == 0xffffffff:
        raise MemoryError

    return indices

def split_labels(indices1, indices2):
    if len(indices1.shape) != 1 or len(indices2.shape) != 1:
        raise ValueError("indices array does not have correct dimensions")
    if indices1.shape[0] != indices2.shape[0]:
        raise ValueError("indices arrays have different length")

    indices1 = np.asarray(indices1.copy(), dtype=np.uint32)
    indices2 = np.asarray(indices2, dtype=np.uint32)
    lib.split_labels(indices1, indices2, indices1.shape[0])
    return indices1

def intersection_matrix(indices1, indices2):
    if len(indices1.shape) != 1 or len(indices2.shape) != 1:
        raise ValueError("indices array does not have correct dimensions")
    if indices1.shape[0] != indices2.shape[0]:
        raise ValueError("indices arrays have different length")

    indices1 = np.asarray(indices1, dtype=np.uint32)
    indices2 = np.asarray(indices2, dtype=np.uint32)
    return Graph(obj=lib.intersection_matrix(indices1, indices2, indices1.shape[0]))

def confusion_matrix(indices_true, indices_pred):
    if len(indices_true.shape) != 1 or len(indices_pred.shape) != 1:
        raise ValueError("indices array does not have correct dimensions")
    if indices_true.shape[0] != indices_pred.shape[0]:
        raise ValueError("indices arrays have different length")

    indices_true = np.asarray(indices_true, dtype=np.uint32)
    indices_pred = np.asarray(indices_pred, dtype=np.uint32)
    a, b, c, d = c_double(), c_double(), c_double(), c_double()
    if not lib.confusion_matrix(a, b, c, d, indices_true, indices_pred, indices_true.shape[0]):
        raise MemoryError

    return np.array([[a.value, b.value], [c.value, d.value]])

def __wrap_metric(fn):
    def wrap(indices_true, indices_pred):
        if len(indices_true.shape) != 1 or len(indices_pred.shape) != 1:
            raise ValueError("indices array does not have correct dimensions")
        if indices_true.shape[0] != indices_pred.shape[0]:
            raise ValueError("indices arrays have different length")

        indices_true = np.asarray(indices_true, dtype=np.uint32)
        indices_pred = np.asarray(indices_pred, dtype=np.uint32)
        return fn(indices_true, indices_pred, indices_true.shape[0])
    return wrap

precision           = __wrap_metric(lib.precision)
recall              = __wrap_metric(lib.recall)
rand_index          = __wrap_metric(lib.rand_index)
fowlkes_mallows     = __wrap_metric(lib.fowlkes_mallows)
jaccard             = __wrap_metric(lib.jaccard)
f1_measure          = __wrap_metric(lib.f1_measure)
adj_rand_index      = __wrap_metric(lib.adj_rand_index)
norm_mutual_info    = __wrap_metric(lib.norm_mutual_info)

if __name__ == '__main__':
    import itertools
    import tempfile
    import unittest
    import random

    adjnoun_maxmodularity = [111, 111, 111, 111, 110, 110, 108, 108, 111, 111, 102, 102, 102, 108,
                             108, 106, 108, 107, 102, 111, 110, 102, 111, 108, 108, 111, 111, 110,
                             107, 74,  107, 107, 107, 107, 107, 107, 108, 108, 108, 110, 111, 111,
                             111, 110, 102, 111, 111, 108, 102, 108, 102, 108, 107, 102, 108, 108,
                             111, 108, 74,  108, 106, 111, 106, 110, 110, 107, 108, 108, 108, 102,
                             108, 106, 108, 111, 74,  108, 106, 111, 102, 110, 110, 110, 107, 110,
                             107, 111, 107, 108, 106, 102, 111, 111, 106, 107, 110, 106, 106, 111,
                             110, 102, 110, 106, 102, 107, 107, 108, 106, 107, 108, 110, 110, 111]

    adjnoun_closeness_centrality =  [0.46835443037974683, 0.37500000000000000, 0.54679802955665020, 0.43700787401574803,
                                     0.42692307692307690, 0.31988472622478387, 0.38275862068965516, 0.44223107569721115,
                                     0.48471615720524020, 0.32743362831858410, 0.34905660377358490, 0.39222614840989400,
                                     0.44578313253012050, 0.44939271255060730, 0.41417910447761190, 0.47033898305084750,
                                     0.34049079754601225, 0.61666666666666670, 0.46058091286307050, 0.39928057553956836,
                                     0.37123745819397990, 0.46835443037974683, 0.41263940520446096, 0.42857142857142855,
                                     0.47844827586206895, 0.47639484978540775, 0.48260869565217390, 0.48898678414096920,
                                     0.44758064516129030, 0.36754966887417220, 0.40959409594095940, 0.47234042553191490,
                                     0.45121951219512196, 0.37755102040816324, 0.44047619047619047, 0.43700787401574803,
                                     0.44578313253012050, 0.42692307692307690, 0.44400000000000000, 0.39084507042253520,
                                     0.35576923076923080, 0.43023255813953487, 0.42528735632183906, 0.51152073732718900,
                                     0.43190661478599224, 0.39084507042253520, 0.36274509803921570, 0.36274509803921570,
                                     0.43359375000000000, 0.39084507042253520, 0.47234042553191490, 0.52606635071090050,
                                     0.42528735632183906, 0.37627118644067800, 0.46835443037974683, 0.33035714285714285,
                                     0.34796238244514105, 0.36393442622950820, 0.32173913043478260, 0.45679012345679010,
                                     0.38144329896907214, 0.35126582278481010, 0.33333333333333330, 0.38013698630136990,
                                     0.27611940298507465, 0.40959409594095940, 0.44400000000000000, 0.41886792452830190,
                                     0.43190661478599224, 0.31444759206798867, 0.44578313253012050, 0.35922330097087380,
                                     0.41729323308270677, 0.40959409594095940, 0.26941747572815533, 0.40363636363636360,
                                     0.41573033707865170, 0.39928057553956836, 0.36038961038961037, 0.40808823529411764,
                                     0.44047619047619047, 0.39501779359430605, 0.38275862068965516, 0.43700787401574803,
                                     0.42528735632183906, 0.39084507042253520, 0.39084507042253520, 0.42528735632183906,
                                     0.41263940520446096, 0.42692307692307690, 0.29919137466307280, 0.31988472622478387,
                                     0.36393442622950820, 0.40217391304347827, 0.32080924855491330, 0.26746987951807230,
                                     0.34365325077399383, 0.35463258785942490, 0.43359375000000000, 0.35806451612903230,
                                     0.35806451612903230, 0.35576923076923080, 0.43529411764705883, 0.44047619047619047,
                                     0.48898678414096920, 0.40808823529411764, 0.38675958188153310, 0.31534090909090910,
                                     0.35576923076923080, 0.29057591623036650, 0.30662983425414364, 0.35463258785942490]

    adjnoun_harmonic_centrality =   [0.52177177177177210, 0.40315315315315325, 0.62912912912912920, 0.48048048048048100,
                                     0.46621621621621640, 0.34009009009008995, 0.42192192192192250, 0.49174174174174200,
                                     0.54204204204204230, 0.34594594594594590, 0.38093093093093110, 0.43753753753753800,
                                     0.49099099099099136, 0.50225225225225260, 0.45195195195195220, 0.51501501501501520,
                                     0.36891891891891910, 0.71021021021021000, 0.51426426426426470, 0.44669669669669700,
                                     0.40915915915915930, 0.50825825825825850, 0.44819819819819834, 0.47072072072072100,
                                     0.53153153153153180, 0.52702702702702730, 0.52177177177177200, 0.53753753753753770,
                                     0.50225225225225260, 0.40315315315315340, 0.44519519519519540, 0.51651651651651680,
                                     0.49474474474474510, 0.41291291291291304, 0.47822822822822847, 0.47222222222222243,
                                     0.48273273273273287, 0.47672672672672695, 0.48723723723723766, 0.42642642642642675,
                                     0.38768768768768797, 0.47147147147147184, 0.46396396396396433, 0.59234234234234250,
                                     0.46996996996997026, 0.43018018018018050, 0.39189189189189190, 0.39789789789789815,
                                     0.48123123123123157, 0.42717717717717750, 0.52702702702702740, 0.60135135135135130,
                                     0.46471471471471500, 0.41366366366366400, 0.51726726726726760, 0.35825825825825830,
                                     0.37642642642642665, 0.39594594594594630, 0.34339339339339325, 0.50975975975976010,
                                     0.41591591591591626, 0.37717717717717730, 0.35720720720720730, 0.41741741741741767,
                                     0.29204204204204180, 0.44819819819819845, 0.48123123123123140, 0.45045045045045060,
                                     0.48573573573573610, 0.33723723723723725, 0.49849849849849860, 0.39894894894894933,
                                     0.45945945945945965, 0.45195195195195226, 0.28483483483483474, 0.45570570570570600,
                                     0.47222222222222240, 0.43018018018018040, 0.38813813813813847, 0.46771771771771780,
                                     0.48048048048048060, 0.43468468468468480, 0.40615615615615625, 0.46771771771771786,
                                     0.46096096096096120, 0.41816816816816830, 0.41666666666666685, 0.47447447447447460,
                                     0.45120120120120144, 0.46696696696696727, 0.31441441441441430, 0.33798798798798796,
                                     0.40240240240240250, 0.42867867867867880, 0.34489489489489494, 0.28363363363363353,
                                     0.37012012012012013, 0.37747747747747745, 0.46546546546546570, 0.38363363363363380,
                                     0.38813813813813836, 0.38093093093093120, 0.48423423423423450, 0.48423423423423450,
                                     0.55705705705705740, 0.45405405405405436, 0.42642642642642680, 0.33573573573573570,
                                     0.38588588588588600, 0.31156156156156145, 0.33168168168168160, 0.37747747747747745]

    class PyGraphTests(unittest.TestCase):
        def test_get_link(self):
            expected = np.zeros(shape=(100,))
            g = Graph(100, directed=True, sorted=True)
            for i in random.sample(range(100), 100):
                expected[i] = 1.0
                g[0, i] = 1.0
                for j in xrange(100):
                    self.assertEqual(g[0, j], expected[j])
            del g

            expected = np.zeros(shape=(100,))
            g = Graph(100, directed=True, sorted=False)
            for i in random.sample(range(100), 100):
                expected[i] = 1.0
                g[0, i] = 1.0
                for j in xrange(100):
                    self.assertEqual(g[0, j], expected[j])
            del g

        def test_load_graph(self):
            g = Graph.load_graph("data/example-trivial.graph")
            self.assertEqual(g.num_nodes, 9)
            self.assertEqual(g.flags, 0)
            self.assertEqual(g.directed, False)
            self.assertEqual(g[0, 0], 0.0)
            self.assertFalse(g.has_edge((0, 0)))
            self.assertEqual(g[0, 1], 1.0)
            self.assertTrue(g.has_edge((0, 1)))
            self.assertEqual(g[0, 2], 0.0)
            self.assertFalse(g.has_edge((0, 2)))
            self.assertEqual(g[1, 0], 1.0)
            self.assertTrue(g.has_edge((1, 0)))
            self.assertEqual(g[1, 1], 0.0)
            self.assertFalse(g.has_edge((1, 1)))
            self.assertEqual(g[1, 2], 1.0)
            self.assertTrue(g.has_edge((1, 2)))
            self.assertEqual(g[2, 0], 0.0)
            self.assertFalse(g.has_edge((2, 0)))
            self.assertEqual(g[2, 1], 1.0)
            self.assertTrue(g.has_edge((2, 1)))
            self.assertEqual(g[2, 2], 0.0)
            self.assertFalse(g.has_edge((2, 2)))
            self.assertEqual(g.count_links(), 13)
            self.assertEqual(g.count_links(directed=True), 26)
            self.assertEqual(g.sum_weights(), 13.0)
            self.assertEqual(g.edges().tolist(), [[0, 1], [0, 3], [0, 4], [1, 2],
                                                  [1, 4], [2, 3], [2, 5], [3, 4],
                                                  [5, 6], [5, 7], [5, 8], [6, 7], [7, 8]])
            g = g.copy()
            self.assertEqual(g.edges().tolist(), [[0, 1], [0, 3], [0, 4], [1, 2],
                                                  [1, 4], [2, 3], [2, 5], [3, 4],
                                                  [5, 6], [5, 7], [5, 8], [6, 7], [7, 8]])
            del g

            with self.assertRaises(IOError):
                Graph.load_graph("does-not-exist.txt")

            with tempfile.NamedTemporaryFile() as fp:
                with self.assertRaises(IOError):
                    Graph.load_graph(fp.name)

            with self.assertRaises(ValueError):
                Graph(obj=0)

        def test_copy(self):
            g = Graph(num_nodes=2, directed=True, sorted=False)
            self.assertEqual(g.flags, 5)
            self.assertEqual(g.directed, True)
            g[0, 1] = 1.0
            g = g.copy()
            self.assertEqual(g[0, 1], 1.0)
            self.assertEqual(g[1, 0], 0.0)
            g.add_edge((0, 1), 4.0)
            self.assertEqual(g[0, 1], 5.0)
            del g[0, 1]
            self.assertEqual(g[0, 1], 0.0)
            del g

        def test_transpose(self):
            g = Graph(num_nodes=2, directed=True)
            g[0, 1] = 1.0
            g = g.transpose()
            self.assertEqual(g[0, 1], 0.0)
            self.assertEqual(g[1, 0], 1.0)
            del g

        def test_invert(self):
            g = Graph(num_nodes=2, directed=True)
            g[0, 1] = 1.0
            h = g.invert()
            self.assertEqual(h[0, 0], 1.0)
            self.assertEqual(h[0, 1], 0.0)
            self.assertEqual(h[1, 0], 1.0)
            self.assertEqual(h[1, 1], 1.0)
            del h
            h = g.invert(self_loops=False)
            self.assertEqual(h[0, 0], 0.0)
            self.assertEqual(h[0, 1], 0.0)
            self.assertEqual(h[1, 0], 1.0)
            self.assertEqual(h[1, 1], 0.0)
            del g
            del h

            g = Graph.load_graph("data/example-trivial.graph")
            h = g.invert()
            self.assertEqual(h.edges().tolist(), [[0, 0], [0, 2], [0, 5], [0, 6], [0, 7], [0, 8],
                                                  [1, 1], [1, 3], [1, 5], [1, 6], [1, 7], [1, 8],
                                                  [2, 2], [2, 4], [2, 6], [2, 7], [2, 8], [3, 3],
                                                  [3, 5], [3, 6], [3, 7], [3, 8], [4, 4], [4, 5],
                                                  [4, 6], [4, 7], [4, 8], [5, 5], [6, 6], [6, 8],
                                                  [7, 7], [8, 8]])
            self.assertEqual(h.weights().tolist(), [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                                    1.0, 1.0, 1.0, 1.0, 1.0])
            del h
            h = g.invert(self_loops=False)
            self.assertEqual(h.edges().tolist(), [[0, 2], [0, 5], [0, 6], [0, 7], [0, 8], [1, 3],
                                                  [1, 5], [1, 6], [1, 7], [1, 8], [2, 4], [2, 6],
                                                  [2, 7], [2, 8], [3, 5], [3, 6], [3, 7], [3, 8],
                                                  [4, 5], [4, 6], [4, 7], [4, 8], [6, 8]])
            self.assertEqual(h.weights().tolist(), [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                                    1.0, 1.0, 1.0, 1.0, 1.0])
            del g
            del h

        def test_multiply_const(self):
            g = Graph(num_nodes=2, directed=True)
            g[0, 0] = 1.0
            g[0, 1] = 2.0
            g[1, 0] = 3.0
            g[1, 1] = 4.0

            g = g.multiply_const(2.0)
            self.assertEqual(g[0, 0], 2.0)
            self.assertEqual(g[0, 1], 4.0)
            self.assertEqual(g[1, 0], 6.0)
            self.assertEqual(g[1, 1], 8.0)
            del g

        def test_multiply_elementwise(self):
            g = Graph(num_nodes=2, directed=True)
            g[0, 0] = 1.0
            g[0, 1] = 2.0
            g[1, 0] = 3.0
            g[1, 1] = 4.0

            h = Graph(num_nodes=2, directed=True)
            h[0, 0] = 5.0
            h[0, 1] = 6.0
            h[1, 0] = 7.0
            h[1, 1] = 8.0

            with self.assertRaises(ValueError):
                g.multiply_elementwise(Graph(num_nodes=3, directed=True))

            g = g.multiply_elementwise(h)
            self.assertEqual(g[0, 0], 5.0)
            self.assertEqual(g[0, 1], 12.0)
            self.assertEqual(g[1, 0], 21.0)
            self.assertEqual(g[1, 1], 32.0)
            del g
            del h

        def test_scalar_product(self):
            g = Graph(num_nodes=2, directed=True)
            g[0, 0] = 1.0
            g[0, 1] = 2.0
            g[1, 0] = 3.0
            g[1, 1] = 4.0

            h = Graph(num_nodes=2, directed=True)
            h[0, 0] = 5.0
            h[0, 1] = 6.0
            h[1, 0] = 7.0
            h[1, 1] = 8.0

            with self.assertRaises(ValueError):
                g.scalar_product(Graph(num_nodes=3, directed=True))

            self.assertEqual(g.scalar_product(h), 5.0 + 12.0 + 21.0 + 32.0)
            del g
            del h

        def test_multiply(self):
            g = Graph(num_nodes=2, directed=True)
            g[0, 0] = 1.0
            g[0, 1] = 2.0
            g[1, 0] = 3.0
            g[1, 1] = 4.0

            h = Graph(num_nodes=2, directed=True)
            h[0, 0] = 5.0
            h[0, 1] = 6.0
            h[1, 0] = 7.0
            h[1, 1] = 8.0

            with self.assertRaises(ValueError):
                g.multiply(Graph(num_nodes=3, directed=True))

            g = g.multiply(h)
            self.assertEqual(g[0, 0], 19.0)
            self.assertEqual(g[0, 1], 22.0)
            self.assertEqual(g[1, 0], 43.0)
            self.assertEqual(g[1, 1], 50.0)
            del g
            del h

        def test_square(self):
            g = Graph(num_nodes=2, directed=True)
            g[0, 0] = 1.0
            g[0, 1] = 2.0
            g[1, 0] = 3.0
            self.assertEqual(g.edges().tolist(), [[0, 0], [0, 1], [1, 0]])
            self.assertEqual(g.weights().tolist(), [1.0, 2.0, 3.0])
            g = g.square()
            self.assertEqual(g[0, 0], 7.0)
            self.assertEqual(g[0, 1], 2.0)
            self.assertEqual(g[1, 0], 3.0)
            self.assertEqual(g[1, 1], 6.0)
            del g

        def test_filter_labels(self):
            g = Graph(num_nodes=2, directed=True)
            g[0, 0] = 1.0
            g[0, 1] = 2.0
            g[1, 0] = 3.0
            g[1, 1] = 4.0
            g = g.filter_labels(np.array([0, 1]))
            self.assertEqual(g[0, 0], 0.0)
            self.assertEqual(g[0, 1], 2.0)
            self.assertEqual(g[1, 0], 3.0)
            self.assertEqual(g[1, 1], 0.0)
            g = g.filter_labels(np.array([0, 0]))
            self.assertEqual(g[0, 0], 0.0)
            self.assertEqual(g[0, 1], 0.0)
            self.assertEqual(g[1, 0], 0.0)
            self.assertEqual(g[1, 1], 0.0)
            del g

        def test_clamp(self):
            g = Graph(num_nodes=2, directed=True)
            g[0, 0] = 1.0
            g[0, 1] = 2.0
            g[1, 0] = 3.0
            g[1, 1] = 4.0
            g = g.clamp(2.0, 3.0)
            self.assertEqual(g[0, 0], 2.0)
            self.assertEqual(g[0, 1], 2.0)
            self.assertEqual(g[1, 0], 3.0)
            self.assertEqual(g[1, 1], 3.0)
            g = g.clamp(1.0, 1.0)
            self.assertEqual(g[0, 0], 1.0)
            self.assertEqual(g[0, 1], 1.0)
            self.assertEqual(g[1, 0], 1.0)
            self.assertEqual(g[1, 1], 1.0)
            del g

        def test_clamp_invalid(self):
            g = Graph(num_nodes=2, directed=True)
            g[0, 0] = 1.0
            g[0, 1] = 2.0
            g[1, 0] = 3.0
            g[1, 1] = 4.0
            g = g.clamp(3.0, 2.0)
            self.assertEqual(g[0, 0], 3.0)
            self.assertEqual(g[0, 1], 3.0)
            self.assertEqual(g[1, 0], 2.0)
            self.assertEqual(g[1, 1], 2.0)
            del g

        def test_resize(self):
            g = Graph(num_nodes=0, directed=True)
            self.assertEqual(g.num_nodes, 0)
            g.resize(num_nodes=1)
            g[0, 0] = 1.0
            g.resize(num_nodes=2)
            g[0, 1] = 2.0
            g[1, 0] = 3.0
            g[1, 1] = 4.0
            self.assertEqual(g.weights().tolist(), [1.0, 2.0, 3.0, 4.0])
            del g

            g = Graph(num_nodes=2, directed=True)
            g[0, 0] = 1.0
            g[0, 1] = 2.0
            g[1, 0] = 4.0
            g[1, 1] = 5.0
            g.resize(num_nodes=3)
            g[0, 2] = 3.0
            g[1, 2] = 6.0
            g[2, 0] = 7.0
            g[2, 1] = 8.0
            g[2, 2] = 9.0
            self.assertEqual(g.weights().tolist(), [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
            g.resize(num_nodes=2)
            self.assertEqual(g.weights().tolist(), [1.0, 2.0, 4.0, 5.0])
            g.resize(num_nodes=1)
            self.assertEqual(g.weights().tolist(), [1.0])
            g.resize(num_nodes=0)
            self.assertEqual(g.num_nodes, 0)
            del g

        def test_modularity(self):
            g = Graph.load_graph("data/example-adjnoun.graph")
            mod = g.modularity(np.array(adjnoun_maxmodularity))
            self.assertAlmostEqual(mod, 0.294696, places=6)
            del g

        def test_modularity_decode(self):
            g = Graph.load_graph("data/example-adjnoun.graph")
            mod = g.modularity_decode(np.array(adjnoun_maxmodularity))
            self.assertAlmostEqual(mod, 0.294696, places=6)
            del g

        def test_antimodularity(self):
            g = Graph.load_graph("data/example-adjnoun.graph")
            mod = g.antimodularity(np.array(adjnoun_maxmodularity))
            self.assertAlmostEqual(mod, 20.298629, places=6)
            del g

        def test_antimodularity_decode(self):
            g = Graph.load_graph("data/example-adjnoun.graph")
            mod = g.antimodularity_decode(np.array(adjnoun_maxmodularity))
            self.assertAlmostEqual(mod, 20.298629, places=6)
            del g

        def test_decode_labels(self):
            for indices in itertools.combinations_with_replacement(range(6), 6):
                indices = np.array(indices)
                labels = decode_labels(indices)

                expected = np.zeros((indices.shape[0],), dtype=int)
                expected[:] = range(indices.shape[0],)
                for i, j in enumerate(indices):
                    expected[np.where(expected == expected[i])] = expected[j]

                d12 = set(zip(expected, labels))
                d1  = set(expected)
                d2  = set(labels)

                self.assertEqual(len(d1), len(d12))
                self.assertEqual(len(d2), len(d12))
                self.assertEqual(count_labels(labels), len(d12))

            labels = decode_labels(np.array([-1, 0, 0]))
            self.assertEqual(labels.tolist(), [0, 0, 0])
            labels = decode_labels(np.array([-1, 0, 1]))
            self.assertEqual(labels.tolist(), [0, 0, 0])
            labels = decode_labels(np.array([-1, 0, 2]))
            self.assertEqual(labels.tolist(), [0, 0, 2])

            labels = decode_labels(np.array([-1, 1, 0]))
            self.assertEqual(labels.tolist(), [0, 1, 0])
            labels = decode_labels(np.array([-1, 1, 1]))
            self.assertEqual(labels.tolist(), [0, 1, 1])
            labels = decode_labels(np.array([-1, 1, 2]))
            self.assertEqual(labels.tolist(), [0, 1, 2])

            labels = decode_labels(np.array([-1, 2, 0]))
            self.assertEqual(labels.tolist(), [0, 0, 0])
            labels = decode_labels(np.array([-1, 2, 1]))
            self.assertEqual(labels.tolist(), [0, 1, 1])
            labels = decode_labels(np.array([-1, 2, 2]))
            self.assertEqual(labels.tolist(), [0, 1, 1])

        def test_simplify_labels(self):
            indices = np.array([5, 5, 7, 7, 5, 7, 8, 8, 8])
            labels = simplify_labels(indices)
            self.assertEqual(labels.tolist(), [0, 0, 1, 1, 0, 1, 2, 2, 2])

        def test_split_labels(self):
            indices1 = np.array([5, 5, 7, 7, 5, 7, 8, 8, 8])
            indices2 = np.array([0, 1, 1, 1, 0, 0, 1, 0, 1])
            labels = split_labels(indices1, indices2)
            self.assertEqual(labels.tolist(), [0, 1, 3, 3, 0, 2, 5, 4, 5])

        # FIXME: Compare with contingency_matrix from sklearn.metrics.cluster.supervised
        def test_intersection_matrix(self):
            indices1 = np.array([5, 5, 7, 7, 5, 7, 8, 8, 8])
            indices2 = np.array([0, 1, 1, 1, 0, 0, 1, 0, 1])
            g = intersection_matrix(indices1, indices2)
            self.assertEqual(g.edges().tolist(), [[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1]])
            self.assertEqual(g.weights().tolist(), [2.0, 1.0, 1.0, 2.0, 1.0, 2.0])
            del g

        def test_confusion_matrix(self):
            indices_true = np.array([1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3, 1, 1, 3, 3, 3])
            indices_pred = np.array([1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3])
            m = confusion_matrix(indices_true, indices_pred)
            self.assertEqual(m[0, 0], 20.0) # TP
            self.assertEqual(m[0, 1], 20.0) # FP
            self.assertEqual(m[1, 0], 24.0) # FN
            self.assertEqual(m[1, 1], 72.0) # TN

        # FIXME: We need more tests to cover all relevant cases
        def test_metrics(self):
            indices_true = np.array([5, 5, 7, 7, 5, 7, 8, 8, 8])
            indices_pred = np.array([0, 1, 1, 1, 0, 0, 1, 0, 1])
            val = precision(indices_true, indices_pred)
            self.assertAlmostEqual(val, 0.1875000, places=6)
            val = recall(indices_true, indices_pred)
            self.assertAlmostEqual(val, 0.3333333, places=6)
            val = rand_index(indices_true, indices_pred)
            self.assertAlmostEqual(val, 0.4722222, places=6)
            val = fowlkes_mallows(indices_true, indices_pred)
            self.assertAlmostEqual(val, 0.2500000, places=6)
            val = jaccard(indices_true, indices_pred)
            self.assertAlmostEqual(val, 0.1363636, places=6)
            val = f1_measure(indices_true, indices_pred)
            self.assertAlmostEqual(val, 0.2400000, places=6)
            val = adj_rand_index(indices_true, indices_pred)
            self.assertAlmostEqual(val, -0.1176470, places=6)
            val = norm_mutual_info(indices_true, indices_pred)
            self.assertAlmostEqual(val, 0.05650554, places=6)

            indices_pred = np.array([5, 5, 7, 7, 5, 7, 8, 8, 8])
            val = precision(indices_true, indices_pred)
            self.assertEqual(val, 1.0)
            val = recall(indices_true, indices_pred)
            self.assertEqual(val, 1.0)
            val = rand_index(indices_true, indices_pred)
            self.assertEqual(val, 1.0)
            val = fowlkes_mallows(indices_true, indices_pred)
            self.assertEqual(val, 1.0)
            val = jaccard(indices_true, indices_pred)
            self.assertEqual(val, 1.0)
            val = f1_measure(indices_true, indices_pred)
            self.assertEqual(val, 1.0)
            val = adj_rand_index(indices_true, indices_pred)
            self.assertEqual(val, 1.0)
            val = norm_mutual_info(indices_true, indices_pred)
            self.assertEqual(val, 1.0)

            for indices_true in (np.array([1, 1, 1, 1, 1, 1, 1, 1, 1]),
                                 np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])):
                for indices_pred in (np.array([1, 1, 1, 1, 1, 1, 1, 1, 1]),
                                     np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])):
                    expected = 1.0 if np.array_equal(indices_true, indices_pred) else 0.0

                    val = precision(indices_true, indices_pred)
                    self.assertEqual(val, expected)
                    val = recall(indices_true, indices_pred)
                    self.assertEqual(val, expected)
                    val = rand_index(indices_true, indices_pred)
                    self.assertEqual(val, expected)
                    val = fowlkes_mallows(indices_true, indices_pred)
                    self.assertEqual(val, expected)
                    val = jaccard(indices_true, indices_pred)
                    self.assertEqual(val, expected)
                    val = f1_measure(indices_true, indices_pred)
                    self.assertEqual(val, expected)
                    val = adj_rand_index(indices_true, indices_pred)
                    self.assertEqual(val, expected)
                    val = norm_mutual_info(indices_true, indices_pred)
                    self.assertEqual(val, expected)

            indices_true = np.array([1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3, 1, 1, 3, 3, 3])
            indices_pred = np.array([1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3])
            val = precision(indices_true, indices_pred)
            self.assertAlmostEqual(val, 0.5000000, places=6)
            val = recall(indices_true, indices_pred)
            self.assertAlmostEqual(val, 0.4545454, places=6)

        def test_bfs(self):
            g = Graph(5, directed=True)
            g[0, 1] = 1.0
            g[1, 2] = 1.0
            g[2, 3] = 1.0
            g[3, 4] = 1.5
            g[2, 4] = 1.5
            self.assertEqual(g.edges().tolist(), [[0, 1], [1, 2], [2, 3], [2, 4], [3, 4]])
            self.assertEqual(g.weights().tolist(), [1.0, 1.0, 1.0, 1.5, 1.5])

            value = g.get_distance_count(100, 0)
            self.assertEqual(value, 0xffffffff)
            value = g.get_distance_weight(100, 0)
            self.assertEqual(value, np.inf)

            value = g.get_distance_count(0, 0)
            self.assertEqual(value, 0)
            value = g.get_distance_weight(0, 0)
            self.assertEqual(value, 0.0)

            value = g.get_distance_count(0, 4)
            self.assertEqual(value, 3)
            value = g.get_distance_weight(0, 4)
            self.assertEqual(value, 3.5)

            counts = g.get_all_distances_count(0)
            self.assertEqual(counts.tolist(), [0, 1, 2, 3, 3])
            counts = g.get_all_distances_count(0, max_count=2)
            self.assertEqual(counts.tolist(), [0, 1, 2, 0xffffffff, 0xffffffff])

            counts = g.get_all_distances_count(100)
            self.assertEqual(counts[0], 0xffffffff)

            weights = g.get_all_distances_weight(0)
            self.assertEqual(weights.tolist(), [0.0, 1.0, 2.0, 3.0, 3.5])
            weights = g.get_all_distances_weight(0, max_weight=2.0)
            self.assertEqual(weights.tolist(), [0.0, 1.0, 2.0, np.inf, np.inf])

            weights = g.get_all_distances_weight(100)
            self.assertEqual(weights[0], np.inf)

            results = g.bfs_count(0, max_count=2)
            self.assertEqual(results, [(0.0, 0, 0xffffffff, 0), (1.0, 1, 0, 1), (2.0, 2, 1, 2)])

            results = g.bfs_count(0)
            self.assertEqual(results, [(0.0, 0, 0xffffffff, 0), (1.0, 1, 0, 1),
                                       (2.0, 2, 1, 2), (3.0, 3, 2, 3), (3.5, 3, 2, 4)])

            results = g.bfs_weight(0, max_weight=2.0)
            self.assertEqual(results, [(0.0, 0, 0xffffffff, 0), (1.0, 1, 0, 1), (2.0, 2, 1, 2)])

            results = g.bfs_weight(0)
            self.assertEqual(results, [(0.0, 0, 0xffffffff, 0), (1.0, 1, 0, 1),
                                       (2.0, 2, 1, 2), (3.0, 3, 2, 3), (3.5, 3, 2, 4)])

            distances = g.get_all_distances_graph(use_weights=False)
            self.assertEqual(distances.edges().tolist(), [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]])
            self.assertEqual(distances.weights().tolist(), [1.0, 2.0, 3.0, 3.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0])
            del distances

            distances = g.get_all_distances_graph(use_weights=True)
            self.assertEqual(distances.edges().tolist(), [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]])
            self.assertEqual(distances.weights().tolist(), [1.0, 2.0, 3.0, 3.5, 1.0, 2.0, 2.5, 1.0, 1.5, 1.5])
            del distances

            del g

        def test_connected_components(self):
            g = Graph.load_graph("data/example-adjnoun.graph")
            components = g.get_connected_components()
            self.assertEqual(len(np.unique(components)), 1)
            del g

        def test_sum_weights(self):
            g = Graph(2, directed=True)
            g[0, 0] = 1.0
            g[0, 1] = 2.0
            g[1, 0] = 3.0

            self.assertEqual(g.sum_weights(), 6.0)

            degrees = g.get_out_degrees()
            self.assertEqual(degrees.tolist(), [2, 1])
            weights = g.get_out_weights()
            self.assertEqual(weights.tolist(), [3.0, 3.0])

            degrees = g.get_in_degrees()
            self.assertEqual(degrees.tolist(), [2, 1])
            weights = g.get_in_weights()
            self.assertEqual(weights.tolist(), [4.0, 2.0])
            del g

        def test_merge_nodes(self):
            g = Graph(4, directed=False)

            c = 0
            for i in xrange(4):
                for j in xrange(i + 1):
                    g[i, j] = c
                    c = c + 1

            g.merge_nodes_custom(0, 2, 0.5, 0.5, 0.0, 0.5)

            # invalid merges
            g.merge_nodes_custom(0, 0, 0.5, 0.5, 0.0, 0.5)
            g.merge_nodes_custom(0, 4, 0.5, 0.5, 0.0, 0.5)
            g.merge_nodes_custom(4, 0, 0.5, 0.5, 0.0, 0.5)

            for i in xrange(4):
                self.assertFalse(g.has_edge((i, 2)))
                self.assertFalse(g.has_edge((2, i)))

            self.assertFalse(g.has_edge((0, 0)))
            self.assertEqual(g[0, 1], 4.0)
            self.assertEqual(g[0, 3], 8.0)
            self.assertEqual(g[1, 0], 4.0)
            self.assertEqual(g[1, 1], 2.0)
            self.assertEqual(g[1, 3], 7.0)
            self.assertEqual(g[3, 0], 8.0)
            self.assertEqual(g[3, 1], 7.0)
            self.assertEqual(g[3, 3], 9.0)
            del g

        def test_newman(self):
            g = Graph.load_graph("data/example-adjnoun.graph")
            labels = g.newman(NEWMAN_FLAGS_FAST | NEWMAN_FLAGS_STRICT, 0)
            self.assertEqual(labels.tolist(), adjnoun_maxmodularity)
            del g

        def test_centrality(self):
            g = Graph.load_graph("data/example-adjnoun.graph")
            centrality = g.closeness_centrality()
            self.assertEqual(centrality.tolist(), adjnoun_closeness_centrality)
            centrality = g.harmonic_centrality()
            self.assertEqual(centrality.tolist(), adjnoun_harmonic_centrality)
            for max_count in [0, 1, 2, 3]:
                centrality_orig = g.harmonic_centrality(max_count=1)
                centrality_fast = g.harmonic_centrality(max_count=1, fast=True)
                self.assertTrue(sum((centrality_orig - centrality_fast)**2) < 1e-10)
            del g

        def test_transitivity(self):
            g = Graph(num_nodes=3, directed=False)
            g[0, 1] = g[0, 2] = 1.0
            g[1, 2] = 1.0
            self.assertEqual(g.transitivity(), 1.0)
            g[0, 0] = 1.0
            self.assertEqual(g.transitivity(), 1.0)
            del g
            g = Graph(num_nodes=5, directed=False)
            g[0, 1] = g[0, 2] = g[0, 3] = g[0, 4] = 1.0
            g[1, 2] = g[1, 3] = g[1, 4] = 1.0
            g[2, 3] = g[2, 4] = 1.0
            g[3, 4] = 1.0
            self.assertEqual(g.transitivity(), 1.0)
            g[0, 0] = 1.0
            self.assertEqual(g.transitivity(), 1.0)
            del g
            g = Graph.load_graph("data/example-adjnoun.graph")
            self.assertEqual(g.transitivity(), 0.15693497881746177)
            del g

        def test_clustering_coefficient(self):
            g = Graph(num_nodes=3, directed=False)
            g[0, 1] = g[0, 2] = 1.0
            g[1, 2] = 1.0
            clustering = g.clustering_coefficient()
            self.assertEqual(clustering.tolist(), [1.0, 1.0, 1.0])
            g[0, 0] = 1.0
            clustering = g.clustering_coefficient()
            self.assertEqual(clustering.tolist(), [1.0, 1.0, 1.0])
            del g
            g = Graph(num_nodes=5, directed=False)
            g[0, 1] = g[0, 2] = g[0, 3] = g[0, 4] = 1.0
            g[1, 2] = g[1, 3] = g[1, 4] = 1.0
            g[2, 3] = g[2, 4] = 1.0
            g[3, 4] = 1.0
            clustering = g.clustering_coefficient()
            self.assertEqual(clustering.tolist(), [1.0, 1.0, 1.0, 1.0, 1.0])
            g[0, 0] = 1.0
            clustering = g.clustering_coefficient()
            self.assertEqual(clustering.tolist(), [1.0, 1.0, 1.0, 1.0, 1.0])
            del g
            g = Graph.load_graph("data/example-adjnoun.graph")
            clustering = g.clustering_coefficient()
            self.assertEqual(sum(clustering), 19.358088938761203)
            self.assertEqual(sum(clustering**2), 6.9377394578170755)
            del g

        def test_reciprocity(self):
            g = Graph(num_nodes=3, directed=True)
            g[0, 1] = g[1, 2] = g[2, 0] = 1.0
            self.assertEqual(g.reciprocity(), 0.0)
            g[1, 0] = 1.0
            self.assertEqual(g.reciprocity(), 2.0 / 4.0)
            g[2, 1] = 1.0
            self.assertEqual(g.reciprocity(), 4.0 / 5.0)
            g[0, 2] = 1.0
            self.assertEqual(g.reciprocity(), 1.0)
            g[0, 0] = 1.0
            self.assertEqual(g.reciprocity(), 1.0)
            del g[0, 2]
            self.assertEqual(g.reciprocity(), 5.0 / 6.0)
            del g

        def test_adj_rand_index_comp(self):
            g = Graph(num_nodes=4, directed=False)
            g[0, 1] = g[2, 3] = 1.0

            indices_true = np.array([0, 1, 0, 1])
            indices_pred = np.array([0, 1, 2, 3])
            val = adj_rand_index(indices_true, indices_pred)
            self.assertEqual(val, 0.0)
            val = g.adj_rand_index_comp(indices_true, indices_pred)
            self.assertEqual(val, 1.0)

            indices_true = np.array([0, 1, 0, 1])
            indices_pred = np.array([0, 1, 0, 3])
            val = adj_rand_index(indices_true, indices_pred)
            self.assertAlmostEqual(val, 0.5714286, places=6)
            val = g.adj_rand_index_comp(indices_true, indices_pred)
            self.assertEqual(val, 1.0)

            indices_true = np.array([0, 1, 0, 1])
            indices_pred = np.array([0, 1, 1, 0])
            val = adj_rand_index(indices_true, indices_pred)
            self.assertEqual(val, -0.5)
            val = g.adj_rand_index_comp(indices_true, indices_pred)
            self.assertEqual(val, 1.0)

            indices_true = np.array([0, 1, 0, 1])
            indices_pred = np.array([0, 1, 1, 1])
            val = adj_rand_index(indices_true, indices_pred)
            self.assertEqual(val, 0.0)
            val = g.adj_rand_index_comp(indices_true, indices_pred)
            self.assertEqual(val, 0.0)
            del g

        def test_batch(self):
            g = Graph(num_nodes=3, directed=True)
            edges = np.array([[0, 1], [1, 2], [2, 0], [3, 3]])
            weights = np.array([1.0, 2.0, 3.0, 4.0])
            g.add_edges(edges, weights)
            self.assertEqual(g.edges().tolist(), [[0, 1], [1, 2], [2, 0]])
            self.assertEqual(g.weights().tolist(), [1.0, 2.0, 3.0])
            g.add_edges(edges, weights)
            self.assertEqual(g.edges().tolist(), [[0, 1], [1, 2], [2, 0]])
            self.assertEqual(g.weights().tolist(), [2.0, 4.0, 6.0])
            g.set_edges(edges, weights)
            self.assertEqual(g.edges().tolist(), [[0, 1], [1, 2], [2, 0]])
            self.assertEqual(g.weights().tolist(), [1.0, 2.0, 3.0])
            g.min_edges(edges, np.array([2.0, 2.0, 2.0, 2.0]))
            self.assertEqual(g.edges().tolist(), [[0, 1], [1, 2], [2, 0]])
            self.assertEqual(g.weights().tolist(), [1.0, 2.0, 2.0])
            g.set_edges(edges, weights)
            g.max_edges(edges, np.array([2.0, 2.0, 2.0, 2.0]))
            self.assertEqual(g.edges().tolist(), [[0, 1], [1, 2], [2, 0]])
            self.assertEqual(g.weights().tolist(), [2.0, 2.0, 3.0])
            g.del_edges(edges)
            self.assertEqual(g.edges().tolist(), [])
            del g

        def test_matrix(self):
            g = Graph(num_nodes=3, directed=True)
            g[0, 2] = 1.0
            g[1, 2] = 2.0
            g[2, 2] = 3.0
            expected = [[0.0, 0.0, 1.0], [0.0, 0.0, 2.0], [0.0, 0.0, 3.0]]
            self.assertTrue(np.array_equal(g.matrix(), expected))
            del g

    # Run the unit tests ...
    unittest.main()
