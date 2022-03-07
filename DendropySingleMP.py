#!/usr/bin/env python3
"""
Bipartition Frequency Hash Robinson Foulds (BFHRF)
In calculating RF between 2 disparate lists of trees, q and r,  we employ a BFH.
This allows us to do a 1 vs BFH(r) for all trees in q.
The normal way is to do i vs j for all i in q and j in r.

"""

import argparse
import dendropy
import numpy as np
from time import time
from multiprocessing import Pool
from math import ceil

__author__ = "Alvin Chon"
__version__ = "0.1.0"
__license__ = "MIT"

# globals
tns = dendropy.TreeList.taxon_namespace
ref_trees_list = []
ref_trees_bipartitions = {}
num_taxa = 0
num_cpu = 0
num_ref_trees = 0


def read_tree(tree):
    """
    Creates Dendropy tree object and encodes bipartitions
    :param tree:
    :return:
    """
    tree_object = dendropy.Tree.get(data=tree, label=tree, schema="newick", taxon_namespace=tns)
    tree_object.encode_bipartitions()
    return tree_object


def parse_trees(tree_files):
    """
    Parses set of newick trees
    :param tree_files:
    :return: tree_list
    """
    tree_list = dendropy.TreeList()
    tree_list.taxon_namespace = tns
    chunk_size = max(ceil(num_cpu / 10), 1)
    pool = Pool(processes=num_cpu)
    for r in pool.imap_unordered(read_tree, tree_files, chunk_size):
        tree_list.append(r)
    pool.close()
    pool.join()
    return tree_list


def rf_single_manual(query_tree_str):
    """
    Computes rf via 1v1 calculation through dendropy of a tree against ref_trees_list
    Manually calculated
    Faster than Dendropy version:  less overhead
    Also, Dendropy's version is not MP friendly.
    :param query_tree_str:
    :return:
    """
    query_tree = dendropy.Tree.get(data=query_tree_str, schema="newick", taxon_namespace=tns)
    query_tree.encode_bipartitions()
    query_bipartition_set = set(query_tree.bipartition_encoding)
    del query_tree
    norm_rf = 0.0
    for ref_trees in ref_trees_list:
        ref_bipartition_set = set(ref_trees.bipartition_encoding)
        rf_left = len(query_bipartition_set.difference(ref_bipartition_set))
        rf_right = len(ref_bipartition_set.difference(query_bipartition_set))
        norm_rf += (rf_left + rf_right)
    return query_tree_str, np.around(norm_rf / num_ref_trees, 5)


def main(args):
    """
    :param args:
    :return:
    """
    global num_taxa
    global num_cpu
    global tns
    global ref_trees_bipartitions
    global ref_trees_list
    global num_ref_trees
    num_taxa = int(args.num_taxa)
    num_cpu = int(args.num_cpu)
    start_time = time()
    begin_time = start_time
    query_trees_rf_dist = {}

    # Sets global tns and gets ref_trees, assumes all trees in q and r have the same tns
    start_time = time()
    ref_trees_files = open(args.ref_trees, 'r').readlines()
    num_ref_trees = len(ref_trees_files)
    print("|Get Files: {}s".format(np.around(time() - start_time, 2)))
    first_tree = dendropy.Tree.get(data=ref_trees_files[0], schema="newick")
    tns = first_tree.taxon_namespace  # sets global taxon_namespace

    # Parse R
    start_time = time()
    ref_trees_list = parse_trees(ref_trees_files)
    num_ref_trees = len(ref_trees_list)
    bipartition_time = time() - start_time
    print('|Parsed {} reference trees and generated bipartition encodings: {}s\tRate:{} trees/s'.format(
        num_ref_trees,
        np.around(bipartition_time, 2),
        np.around(num_ref_trees / bipartition_time, 2)))

    # MP Single list vs list - manually
    # print('|RF single')
    start_time = time()
    trees_rf_single_manual = {}
    chunk_size = max(ceil(num_cpu / 10), 1)
    pool = Pool(processes=num_cpu)
    counter = 0
    query_trees_files = open(args.query_trees, 'r').readlines()
    orig_num_trees = len(query_trees_files)
    total_comp = 1e6
    if len(query_trees_files) * len(ref_trees_list) > total_comp:
        print('|Reducing data size of query list and using estimated rf time for speedup comparison.')
        subset = int(total_comp / len(ref_trees_list))
        if subset <= len(query_trees_files):
            query_trees_list = query_trees_files[0:subset]
        else:
            query_trees_list = query_trees_files[0:1000]
    else:
        query_trees_list = query_trees_files
    for r in pool.imap_unordered(rf_single_manual, query_trees_list, chunk_size):
        counter += 1
        if counter % int(len(query_trees_list) / 10) == 0:
            cur_time = time() - start_time
            trees_per_min = np.around(counter / (cur_time / 60), 2)
            total_rf_time = np.around(orig_num_trees / trees_per_min, 2)
            print("||RF single: {}/{} vs {} in {}s\tAverage(total): {} trees/m\tEstimated total RF time: {}m".format(
                counter, len(query_trees_list), num_ref_trees, np.around(cur_time, 3),
                trees_per_min, total_rf_time))
        tree_name, avg_norm_rf = r
        trees_rf_single_manual[tree_name] = avg_norm_rf
    pool.close()
    pool.join()
    single_time = time() - start_time
    print(
        '|RF single (list vs list): Computed RF of {} query_trees against {} ref_trees: {}s'.format(
            len(query_trees_list),
            num_ref_trees,
            np.around(
                single_time,
                2)))

    # write file
    start_time = time()
    ofh = open(args.output_file, 'w')
    for key in sorted(query_trees_rf_dist.keys()):
        ofh.write("{},{}\n".format(key, str(query_trees_rf_dist[key])))
    ofh.close()
    file_time = time() - start_time
    print('|File Output: {}s'.format(np.around(file_time, 2)))
    #
    print('|STATS: TOTAL TIME: {}s'.format(np.around(time() - begin_time, 2)))
    print('|STATS: TOTAL ESTIMATED TIME: {}m'.format(np.around((bipartition_time + file_time) / 60 + total_rf_time, 2)))


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("ref_trees",
                        help="Required gene tree file.  Assumes one newick tree per line.")
    parser.add_argument("query_trees",
                        help="Required gene tree file.  Assumes one newick tree per line.")
    parser.add_argument("num_taxa", help="Number of taxa")
    parser.add_argument("num_cpu", help="Number of cpus")
    parser.add_argument("-output_file", help="Output file, default=output.txt", default="output.txt")

    input_args = parser.parse_args()
    print(input_args)
    main(input_args)
