# BFHRF
> A command line utility to compute the most parsimonious tree(s) given a collection of Query (Q) and Reference (R) trees for a given number of species (n).  The output is the average traditional Robinson-Foulds (RF) value for each query tree against the reference collection.

## Table of Contents

* [Features](#features)
* [Installation](#installation)
* [Usage](#usage)
* [Data](#data) 
* [Reference](#reference)
* [License](#license)
* [Issues](#issues)

## Features
- Computes average RF value for each query tree against the reference trees.
- Scalable and Extensible:  O(nq + nr) in time and O(nr) in memory
- Uses Python 3.5+ with Dendropy and Multiprocessing packages.
- Easy parallelization with the number of CPUs specified at runtime.

## Installation
Dependency requirements:
- Python 2.6+ or 3.6+
- Dendropy 4+
- multiprocessing

## Usage
BFHRF is intended as a command line utility and not a python package to be easily used in other pipelines.  
Download or clone and run the self-contained python file.
```
$ python3 BipartitionFrequencyHash.py -h
usage: BipartitionFrequencyHash.py [-h] [-output_file OUTPUT_FILE] [-bipartition_filter BIPARTITION_FILTER] ref_trees query_trees num_taxa num_cpu

Visit https://www.github.com/achon/bfhrf/ for more information

positional arguments:
  ref_trees             Required reference tree file. Assumes one newick tree per line.
  query_trees           Required query tree file. Assumes one newick tree per line.
  num_taxa              Number of taxa
  num_cpu               Number of CPUs

optional arguments:
  -h, --help            show this help message and exit
  -output_file OUTPUT_FILE
                        Output file, default=output.txt
  -bipartition_filter BIPARTITION_FILTER
                        Optional bipartition filtering by size. Value can be 1 to floor(n/2). Values entered as range: min-max and does not include max

```
Average RF distance of each tree vs all trees (newick) in a file, trees.tre, all with the same 10 taxa.  The job below uses 3 threads.

`$ python3 BipartitionFrequencyHash.py trees.tre trees.tre 10 3`

## Data
Included are all the data files used in the paper.  See the paper for details.
> [Data and scripts from the paper](https://github.com/achon/bfhrf_data)


## Reference
> [Scalable and Extensible Robinson-Foulds for Comparative Phylogenetics, Chon et al, accepted HiCOMB 2022]

## License
You are free to use the software in any way.  
> [MIT license](http://opensource.org/licenses/mit-license.php)

## Issues
To report bugs, visit the [Issues](https://github.com/achon/bfhrf/issues) page.