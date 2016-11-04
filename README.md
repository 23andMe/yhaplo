# yHaplo: Identifying Y-Chromosome Haplogroups

David Poznik  
23andMe  
October, 2016

--------------------------------------------------------------------------------
## Overview

yHaplo identifies the Y-chromosome haplogroup of each male in a sample of 
one to millions. It does not rely on any particular genotyping modality or platform. 
Although full sequences yield the most granular haplogroup classifications, genotyping 
arrays can yield reliable calls, provided a reasonable number of phylogenetically 
informative variants has been assayed. 

Briefly, haplogroup calling involves two steps. The program first builds an internal 
representation of the Y-chromosome phylogeny by reading its primary structure from 
(Newick-formatted) text and then importing phylogenetically informative SNPs from the 
[ISOGG database](http://isogg.org/tree/ISOGG_YDNA_SNP_Index.html), storing each SNP
within a specific node and growing the tree as necessary. It then traverses the tree
for each individual, identifying for each the path of derived alleles leading to 
a haplogroup designation.

yHaplo is available for non-commercial use pursuant to the terms of the non-exclusive 
license agreement, `LICENSE.txt`. To learn more about the algorithm, please see the 
[white paper](https://api.23andme.com/res/pdf/23-13_paternal_haplogroups_yHaplo.pdf), 
and to learn more about the software, please see `yHaplo.manual.pdf`. 


--------------------------------------------------------------------------------
## Input

### Phylogenetic data

`input/`

* `y.tree.primary.DATE.nwk`   : primary structure of the Y-chromosome tree
* `isogg.DATE.txt`            : phylogenetically informative SNPs
* `isogg.correct.*.txt`       : corrections to ISOGG data
* `isogg.omit.*.txt`          : SNPs to drop due to inconsistencies observed in test data
* `isogg.multiallelic.txt`    : physical coordinates of multiallelic sites to be excluded
* `representative.SNPs.*.txt` : SNPs deemed representative of a corresponding haplogroup


### Supported genotype formats

* `.genos.txt`    : sample-major genotypes  
    * row 1: physical coordinates  
    * column 1: individual IDs
    * cell (i, j): genotype for individual i at position j, encoded as a single character from the set { A, C, G, T, . }, with "." representing an unobserved value
* `.resid.txt`    : file with 23andMe research IDs in the first column
* `.vcf`, `.vcf.gz` : snp-major VCF file
* `.vcf4`         : snp-major pseudo-VCF. differences include:
    * no "#" in header row
    * fewer header columns
    * GT values recorded as { A, C, G, T, . } rather than { 0, 1, . }


--------------------------------------------------------------------------------
## Code

### Driver script

`callHaplogroups.py` : for an overiew of command-line options, issue the following command: `callHaplogroups.py -h`

### Main classes

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
