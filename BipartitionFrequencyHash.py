#!/usr/bin/env python3
"""
Bipartition Frequency Hash Robinson Foulds (BFHRF)
In calculating RF between 2 disparate lists of trees, Q and R,  we employ a BFH.
This allows us to do a q vs BFH_R for all trees in Q.
"""

import argparse
import dendropy
import numpy as np
from time import time
from multiprocessing import Pool
from math import ceil

__author__ = "Alvin Chon"
__email__ = "achon@iastate.edu"
__version__ = "0.1.0"
__license__ = "MIT"

# globals
tns = dendropy.TreeList.taxon_namespace
ref_trees_bipartitions = {}
num_taxa = 0
num_cpu = 0
num_ref_trees = 0
ref_trees_keys = []
ref_trees_sum = 0
bipartition_range = []


# def parse_tree_size(tree):
#     """
#     Parse a tree (Newick str), generate bipartitions, filter based upon sizes, and return bipartitions.
#     :param tree:
#     :return:
#     """
#     tree_object = dendropy.Tree.get(data=tree, schema="newick", taxon_namespace=tns)
#     tree_object.encode_bipartitions()
#     bps = []
#     size_check = True
#     if filter_list == [i for i in range(2, ceil(num_taxa / 2) + 1)]:
#         size_check = False
#     for bbp in tree_object.encode_bipartitions():
#         bp = bin(bbp.split_bitmask)[2:].zfill(num_taxa)
#         if size_check:
#             if bp.count('0') <= bp.count('1'):
#                 char_count = bp.count('0')
#             else:
#                 char_count = bp.count('1')
#             if char_count in filter_list:
#                 bps.append(bp)
#         else:
#             bps.append(bp)
#     return bps

def parse_tree(tree):
    """
    Parse a tree (Newick str) and return Dendropy tree object w/ bipartitions.
    :param tree:
    :return:
    """
    tree_object = dendropy.Tree.get(data=tree, schema="newick", taxon_namespace=tns)
    tree_object.encode_bipartitions()
    return tree_object.bipartition_encoding


def create_bipartition_set(ref_trees_files):
    """
    Parse trees dynamically and build the bfh.
    :param ref_trees_files:
    :return:
    """
    global ref_trees_bipartitions
    chunk_size = max(ceil(num_cpu / 10), 1)
    pool = Pool(processes=num_cpu)
    for r in pool.imap_unordered(parse_tree, ref_trees_files, chunk_size):
        for bp in r:
            try:
                ref_trees_bipartitions[bin(bp.split_bitmask)[2:].zfill(num_taxa)] += 1
            except KeyError:
                ref_trees_bipartitions[bin(bp.split_bitmask)[2:].zfill(num_taxa)] = 1
    pool.close()
    pool.join()
    return


def rf_bfh(query_trees_files):
    """
    Computes RF of Query trees dynamically against BFH_R.
    :param query_trees_files:
    :return:
    """
    global ref_trees_keys
    global ref_trees_sum
    tree_rf_dist = {}
    ref_trees_keys = ref_trees_bipartitions.keys()
    ref_trees_sum = sum(ref_trees_bipartitions.values())
    start_time = time()
    chunk_size = max(ceil(num_cpu / 10), 1)
    pool = Pool(processes=num_cpu)
    for r in pool.imap_unordered(rf_bfh_mp, query_trees_files, chunk_size):
        tree_name, avg_norm_rf = r
        tree_rf_dist[tree_name] = avg_norm_rf
    pool.close()
    pool.join()
    end_time = time() - start_time
    trees_per_min = np.around(len(query_trees_files) / (end_time / 60), 2)
    total_rf_time = np.around(len(query_trees_files) / trees_per_min, 2)
    print("||BHRF: Average(total): {} trees/m\tEstimated total RF time: {}m".format(
        trees_per_min, total_rf_time))
    return tree_rf_dist


def rf_bfh_mp(tree_str):
    """
    Body for computing RF of a tree against ref_trees BFH using MP
    Done using tree1 keys since there are n-1 keys whereas bfh has a minimum of n-1
    bipartition size_filtering included
    :param tree_str:
    :return:
    """
    tree1 = dendropy.Tree.get(data=tree_str, label=tree_str, schema="newick", taxon_namespace=tns)
    tree1.encode_bipartitions()
    tree1_bp = {}
    for bp in tree1.bipartition_encoding:
        bitmask = bin(bp.split_bitmask)[2:].zfill(num_taxa)
        if len(bipartition_range) > 0:
            if bipartition_range[0] <= sum([int(c) for c in bitmask]) <= bipartition_range[1]:
                tree1_bp[bitmask] = num_ref_trees
        else:
            tree1_bp[bitmask] = num_ref_trees
    tree1_keys = tree1_bp.keys()
    # # Sym Diff Left
    # rf_left = sum(tree1_bp.values())
    # for key1 in tree1_keys:
    #     if key1 in ref_trees_keys:
    #         rf_left -= ref_trees_bipartitions[key1]  # tree1 is num_ref_trees where tree2 is at most num_ref_trees
    # # Sym Diff Right
    # rf_right = ref_trees_sum
    # for key1 in tree1_keys:
    #     if key1 in ref_trees_keys:
    #         rf_right -= tree1_bp[key1]
    # Joined Sym Diff
    rf_left = sum(tree1_bp.values())
    rf_right = ref_trees_sum
    for key1 in tree1_keys:
        if key1 in ref_trees_keys:
            rf_left -= ref_trees_bipartitions[key1]  # tree1 is num_ref_trees where tree2 is at most num_ref_trees
            rf_right -= tree1_bp[key1]
    # note that the below is NOT divided by 2
    return tree1.label, (rf_left + rf_right) / num_ref_trees


def main(args):
    """
    :param args:
    :return:
    """
    global num_taxa
    global num_cpu
    global tns
    global ref_trees_bipartitions
    global num_ref_trees
    global bipartition_range
    num_taxa = int(args.num_taxa)
    num_cpu = int(args.num_cpu)
    begin_time = time()

    # Sets global tns and gets ref_trees, assumes all trees in q and r have the same tns
    start_time = begin_time
    ref_trees_files = open(args.ref_trees, 'r').readlines()
    num_ref_trees = len(ref_trees_files)
    print("|Get Reference Trees Files: {}s".format(np.around(time() - start_time, 2)))
    first_tree = dendropy.Tree.get(data=ref_trees_files[0], schema="newick")
    tns = first_tree.taxon_namespace  # sets global taxon_namespace

    # Dynamically read and fill BFH
    create_bipartition_set(ref_trees_files)
    bipartition_time = time() - start_time
    print('|Parsed {} reference trees, generated bipartitions, and created bfh: {}s\tRate:{} trees/s'.format(
        num_ref_trees,
        np.around(bipartition_time, 2),
        np.around(num_ref_trees / bipartition_time, 2)))
    del ref_trees_files  # explicit memory save

    # Get query trees.
    start_time = time()
    query_trees_files = open(args.query_trees, 'r').readlines()
    print("|Get Query Trees Files: {}s".format(np.around(time() - start_time, 2)))

    # Parse, bipartitions, RF vs BFH calc.
    start_time = time()
    query_trees_rf_dist = rf_bfh(query_trees_files)
    bfh_time = time() - start_time
    print(
        '|BFHRF: Computed RF of {} query_trees against {} ref_trees: {}s'.format(len(query_trees_rf_dist.keys()),
                                                                                 num_ref_trees,
                                                                                 np.around(bfh_time, 2)))
    print('|Number of unique bipartitions: {}'.format(len(ref_trees_keys)))
    # write file
    start_time = time()
    ofh = open(args.output_file, 'w')
    for key in sorted(query_trees_rf_dist.keys()):
        ofh.write("{},{}\n".format(key, str(query_trees_rf_dist[key])))
    ofh.close()
    file_time = time() - start_time
    print('|File Output: {}s'.format(np.around(file_time, 2)))
    print('|STATS: TOTAL TIME: {}s'.format(np.around(time() - begin_time, 2)))


if __name__ == "__main__":
    """Bipartition Frequency Hash Robinson Foulds (BFHRF).  In calculating RF between 2 disparate lists of trees, 
    Q and R, we employ a BFH that computes a q vs BFH_R for all trees in Q.  Bipartition filtering is implemented 
    here,  however it has been moved to another repo. It will be available in BPSTERF. """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("ref_trees",
                        help="Required reference tree file.  Assumes one newick tree per line.")
    parser.add_argument("query_trees",
                        help="Required query tree file.  Assumes one newick tree per line.")
    parser.add_argument("num_taxa", help="Number of taxa")
    parser.add_argument("num_cpu", help="Number of CPUs")
    parser.add_argument("-output_file", help="Output file, default=bfhrf_output.txt", default="bfhrf_output.txt")
    # parser.add_argument("-bipartition_filter", help="Optional bipartition filtering by size.  Value can be 1 to "
    #                                                 "floor(n/2).  Values entered as range: min-max and does not "
    #                                                 "include max", type=str, default=[])
    input_args = parser.parse_args()
    print(input_args)
    main(input_args)
