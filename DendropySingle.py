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

__author__ = "Alvin Chon"
__version__ = "0.1.0"
__license__ = "MIT"

# globals
tns = dendropy.TreeList.taxon_namespace
num_taxa = 0
num_ref_trees = 0
orig_len_query_trees = 0


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
    for tree in tree_files:
        tree_list.append(read_tree(tree))
    return tree_list


def rf_single_manual(query_tree, ref_tree):
    """
    Computes rf via 1v1 calculation through dendropy of a tree against ref_trees_list
    Manually calculated
    Faster than Dendropy version:  less overhead
    Also, Dendropy's version is not MP friendly.
    :param query_tree:
    :param ref_tree: 
    :return:
    """
    query_bipartition_set = set(query_tree.bipartition_encoding)
    ref_bipartition_set = set(ref_tree.bipartition_encoding)
    rf_left = len(query_bipartition_set.difference(ref_bipartition_set))
    rf_right = len(ref_bipartition_set.difference(query_bipartition_set))
    norm_rf = (rf_left + rf_right) / (2 * num_taxa - 6)
    return norm_rf


def rf_comparison(trees_list_files, trees_list2):
    """
    Computes avg_norm_rf for all comparisons
    :param trees_list_files:
    :param trees_list2:
    :return:
    """
    start_time = time()
    query_trees_rf_dist = {}
    counter = 0
    for tree_file in trees_list_files:
        tree1 = dendropy.Tree.get(data=tree_file, label=tree_file, schema="newick", taxon_namespace=tns)
        tree1.encode_bipartitions()
        avg_norm_rf = 0
        for tree2 in trees_list2:
            avg_norm_rf += rf_single_manual(tree1, tree2)
        query_trees_rf_dist[tree1.label] = avg_norm_rf
        counter += 1
        if counter % 100 == 0:
            cur_time = time() - start_time
            trees_per_min = np.around(counter / (cur_time / 60), 2)
            total_rf_time = np.around(orig_len_query_trees / trees_per_min, 2)
            print("||RF single: {}/{} vs {} in {}s\tAverage(total): {} trees/m\tEstimated total parse & RF time: {}m".format(
                counter, len(trees_list_files), num_ref_trees, np.around(cur_time, 2),
                trees_per_min, total_rf_time))
            # get out of loop if taking too long
            if cur_time / 60 > 5:
                print("Breaking out of loop...")
                break
    return query_trees_rf_dist, total_rf_time


def main(args):
    """
    :param args:
    :return:
    """
    global num_taxa
    global tns
    global num_ref_trees
    global orig_len_query_trees
    num_taxa = int(args.num_taxa)
    start_time = time()
    begin_time = start_time

    # Sets global tns and gets ref_trees, assumes all trees in q and r have the same tns
    start_time = time()
    ref_trees_files = open(args.ref_trees, 'r').readlines()
    print("|Get Files: {}s".format(np.around(time() - start_time, 2)))
    first_tree = dendropy.Tree.get(data=ref_trees_files[0], schema="newick")
    tns = first_tree.taxon_namespace  # sets global taxon_namespace
    ref_trees_list = parse_trees(ref_trees_files)
    num_ref_trees = len(ref_trees_list)
    del ref_trees_files
    bipartition_time = time() - start_time
    print('|Parsed {} reference trees, bipartition encodings: {}s\tRate:{} trees/s\t'.format(
        num_ref_trees,
        np.around(bipartition_time, 2),
        np.around(num_ref_trees / bipartition_time, 2)))

    # Compute RF of query_trees against ref_trees
    query_trees_files = open(args.query_trees, 'r').readlines()
    orig_len_query_trees = len(query_trees_files)
    if len(query_trees_files) > 1000:
        print("|Trimming list to 1000 files and doing estimation.")
        query_trees_files = query_trees_files[0:1000]
    # query_trees_list = parse_trees(query_trees_files)
    # del query_trees_files
    query_trees_rf_dist, total_rf_time = rf_comparison(query_trees_files, ref_trees_list)

    # write file
    start_time = time()
    ofh = open(args.output_file, 'w')
    for key in sorted(query_trees_rf_dist.keys()):
        ofh.write("{},{}\n".format(key, str(query_trees_rf_dist[key])))
    ofh.close()
    file_time = time() - start_time
    print('|File Output: {}s'.format(np.around(file_time, 2)))
    print('|STATS: TOTAL TIME: {}s'.format(np.around(time() - begin_time, 3)))
    print('|STATS: TOTAL ESTIMATED TIME: {}m'.format(np.around((bipartition_time + file_time) / 60 + total_rf_time, 2)))


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("ref_trees",
                        help="Required gene tree file.  Assumes one newick tree per line.")
    parser.add_argument("query_trees",
                        help="Required gene tree file.  Assumes one newick tree per line.")
    parser.add_argument("num_taxa", help="Number of taxa")
    parser.add_argument("-output_file", help="Output file, default=output.txt", default="output.txt")
    input_args = parser.parse_args()
    print(input_args)
    main(input_args)
