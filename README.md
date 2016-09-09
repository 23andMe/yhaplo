# Y-Chromosome Haplogroup Calling

David Poznik  
23andMe
May, 2016

--------------------------------------------------------------------------------
## Software overview

This software identifies the Y-chromosome haplogroup of each male in a sample of 
one to millions. Sequence data will yield the most highly resolved classifications, 
but the algorithm also works well with array-based genotype data, provided a reasonable 
number of phylogenetically informative sites have been assayed. 

The 2-step process is to:

1. Build internal representation of the Y-chromosome phylogeny.
   a. Read primary structure from text file (Newick format).
   b. Import phylogenetically informative SNPs from the 
      [ISOGG database](http://isogg.org/tree/ISOGG_YDNA_SNP_Index.html),
      storing them within the nodes and growing the tree as necessary.
2. Call haplogroup of each sample.
   a. Build genotype dictionary.
   b. Traverse tree to identify path of derived alleles 
      leading to haplogroup designation.


--------------------------------------------------------------------------------
## Input and options

### Phylogenetic data

input/

* y.tree.primary.<DATE>.nwk  : primary structure of the Y-chromosome tree
* isogg.<DATE>.txt           : phylogenetically informative SNPs
* isogg.correct.*.txt        : corrections to ISOGG data
* isogg.omit.*.txt           : SNPs to drop due to inconsistencies observed in test data
* isogg.multiallelic.txt     : physical coordinates of multiallelic sites to be excluded


### Supported genotype formats

* .genos.txt    : sample-major genotypes  
    * row 1: physical coordinates  
    * column 1: sample IDs
    * cell (i, j): genotype for individual i at position j, encoded as a single character from the set { A, C, G, T, . }, with "." representing an unobserved value
* .resid.txt    : file with 23andMe research IDs in the first column
* .vcf, .vcf.gz : snp-major VCF file
* .vcf4         : snp-major pseudo-VCF. differences include:
    * no "#" in header row
    * fewer header columns
    * GT values recorded as { A, C, G, T, . } rather than { 0, 1, . }

### Command-line options

Please issue the following command to view the full list of command-line options: `callHaplogroups.py -h`


--------------------------------------------------------------------------------
## Code overview

### Main engine

`callHaplogroups.py` : builds tree and calls haplogroups for individuals. 
                          this script's numerous options are described below.

### Classes

* `Tree`         : knows root, depth, haplogroup-to-node mappings, etc.;
                     parses a Newick file to build primary tree;
                     parses ISOGG table to add SNPs to nodes and grow tree;
                     finds the derived path leading from the root to an individual
* `Node`         : element of the tree. knows parent, children, snps, etc.
                    represents the branch that leads to it
* `SNP`          : knows position, ancestral and derived alleles, node, etc.
* `PlatformSNP` : knows position and ablock index 
* `Sample`       : knows an individual's genotypes and haplogroup
* `Customer`     : (subclass of Sample) has 23andMe metadata and genotypes from ablocks
* `Path`         : path through a tree; stores the next node to visit, a list of SNPs 
                    observed in the derived state, the most derived SNP observed, 
                    and the number of ancestral alleles encountered
* `Page`         : 23andMe content page labels
* `Config`       : container for parameters, command-line options, and filenames

### Utilities

* `utils.py`    : shared utility functions

### Auxiliary scripts

* `convert2genos.py` : converts data to .genos.txt format
* `plotTree.py`       : plots a newick tree


--------------------------------------------------------------------------------
## Output directories

* converted/   : output of convert2genos.py
* output/      : output of callHaplogroups.py
