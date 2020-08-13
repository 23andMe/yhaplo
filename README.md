# `yhaplo` | Identifying Y-Chromosome Haplogroups

David Poznik  
23andMe  
October, 2016

--------------------------------------------------------------------------------
## Overview

`yhaplo` identifies the Y-chromosome haplogroup of each male in a sample of one to 
millions. It does not rely on any particular genotyping modality or platform, and it is 
robust to missing data, genotype errors, mutation recurrence, and other complications. 
Although full sequences yield the most granular haplogroup classifications, genotyping 
arrays can yield reliable calls, provided a reasonable number of phylogenetically 
informative variants has been assayed. 

Briefly, haplogroup calling involves two steps. The program first builds an internal 
representation of the Y-chromosome phylogeny by reading its primary structure from 
(Newick-formatted) text and importing phylogenetically informative SNPs from the 
[ISOGG database](http://isogg.org/tree/ISOGG_YDNA_SNP_Index.html), affiliating each 
SNP with the appropriate node and growing the tree as necessary. It then traverses the 
tree for each individual, identifying the path of derived alleles leading to a
haplogroup designation.

`yhaplo` is available for non-commercial use pursuant to the terms of the non-exclusive 
license agreement, `LICENSE.txt`. To learn more about the algorithm, please see our 
bioRxiv [pre-print](http://biorxiv.org/content/early/2016/11/19/088716):

    Poznik GD. 2016. Identifying Y-chromosome haplogroups in arbitrarily large samples 
    of sequenced or genotyped men. bioRxiv doi: 10.1101/088716

To learn more about the software, please see the manual, `yhaplo.manual.<DATE>.pdf`. 
For an overiew of command-line options, install the package and run: `yhaplo -h`


--------------------------------------------------------------------------------
## Installation

`yhaplo` is compatible with both Python 2 and Python 3.  

To install:

```
git clone git@github.com:23andMe/yhaplo.git
cd yhaplo
pip install --editable .
```

To test-run on example data:

```
yhaplo -ex
```

The `-ex` option tells `yhaplo` to run on a subset of 1000 Genomes data
and sets the `--all_aux_output` flag to produce all auxiliary output.


--------------------------------------------------------------------------------
## Caveats

Please note the following caveats before running `yhaplo`:

* `yhaplo` does not check for sex status; it assumes all individuals are male.
* `yhaplo` expects SNP coordinates consistent with the hg19/GRCh37 reference assembly.
* `yhaplo` expects data at a reasonable number of ISOGG SNPs. This assumption is violated by:
  * variants-only sequence data
  * very low-coverage sequencing
  * genotyping arrays with few Y-chromosome probes


If, for a given individual, `yhaplo` observes no derived alleles at ISOGG SNPs on the upper 
branches of the Y-chromosome phylogeny, it will call the individual haplogroup "A," 
since all human Y-chromosome lineages are technically sublineages of A. 
Before concluding that the individual sample belongs to paragroup A (which 
includes haplogroups A00, A0, A1a, and A1b1), run with the `-as` option, and check the 
auxiliary output for ancestral alleles at haplogroup-BT SNPs. If you do not see any, 
your data set probably violates one or more of the assumptions listed above.

In particular, "variants-only" VCF files restrict to SNPs at which alternative alleles 
were observed, but ref/alt status is unimportant to `yhaplo`. What is important is 
ancestral/derived status. The reference sequence contains many derived alleles, 
and `yhaplo` will not be happy if you discard these valuable data. So please emit all 
confident sites when calling variants. To limit compute time and file size, you could 
safely restrict to positions in `output/isogg.snps.unique.DATE.txt`, as these are the 
only SNPs `yhaplo` considers. To generate this file, just run `yhaplo` with no arguments.


--------------------------------------------------------------------------------
## Input

### Phylogenetic data

`input/`

* `y.tree.primary.DATE.nwk` : primary structure of the Y-chromosome tree
* `isogg.DATE.txt` : phylogenetically informative SNPs scraped directly from the ISOGG website. 
`yhaplo` resolves errors and formatting inconsistencies and emits cleaned versions 
(`output/isogg.snps.cleaned.DATE.txt` and `output/isogg.snps.unique.DATE.txt`; 
see `yhaplo.manual.<DATE>.pdf` for details).
* `isogg.correct.*.txt` : corrections to ISOGG data
* `isogg.omit.*.txt` : SNPs to drop due to inconsistencies observed in test data
* `isogg.multiallelic.txt` : physical coordinates of multiallelic sites to be excluded
* `representative.SNPs.*.txt` : SNPs deemed representative of corresponding haplogroups


### Supported genotype formats

* `.genos.txt` : sample-major genotypes  
    * row 1: physical coordinates  
    * column 1: individual IDs
    * cell (i, j): genotype for individual i at position j, encoded as a single character 
from the set { A, C, G, T, . }, with "." representing an unobserved value
* `.resid.txt` : file with 23andMe research IDs in the first column
* `.vcf`, `.vcf.gz` : snp-major VCF file.
    It is most efficient to restrict input VCF files to the Y chromosome.
* `.vcf4` : snp-major pseudo-VCF. differences include:
    * no "#" in header row
    * fewer header columns
    * GT values recorded as { A, C, G, T, . } rather than { 0, 1, . }


--------------------------------------------------------------------------------
## Output

All output file formats are described in detail in `yhaplo.manual.<DATE>.pdf`.

The two primary output files are:

1. `log.projectName.txt` : log file containing details of the run
2. `haplogroups.projectName.txt` : haplogroup calls. The 4 columns are:
    1. ID
    2. Haplogroup short form, with the name of a SNP observed in the derived state
    3. Haplogroup short form, with the name of a representative SNP
    4. Haplogroup long form, using Y-Chromosome Consortium nomenclature

`yhaplo` also produces a number of SNP tables, tree files, and auxiliary output files. 
Please see `yhaplo.manual.<DATE>.pdf` and `yhaplo -h` for details.


--------------------------------------------------------------------------------
## Code

### Driver script

`call_haplogroups.py`

### Main classes

* `Tree` : knows root, depth, haplogroup-to-node mappings, etc.;
          parses a Newick file to build primary tree;
          parses ISOGG table to add SNPs to nodes and grow tree;
          finds the derived path leading from the root to an individual
* `Node` : element of the tree. knows parent, children, snps, etc.
          represents the branch that leads to it
* `SNP` : knows position, ancestral and derived alleles, node, etc.
* `PlatformSNP` : knows position and ablock index 
* `Sample` : knows genotypes and haplogroup of an individual 
* `Customer` : (subclass of Sample) has 23andMe metadata and genotypes from ablocks
* `Path` : path through a tree; stores the next node to visit, a list of SNPs 
          observed in the derived state, the most derived SNP observed, 
          and the number of ancestral alleles encountered
* `Page` : 23andMe content page labels
* `Config` : container for parameters, command-line options, and filenames

### Utilities

* `utils.py` : shared utility functions

### Auxiliary scripts

* `convert_to_genos.py` : converts data to `.genos.txt` format
* `plot_tree.py` : plots a newick tree
