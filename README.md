# Yhaplo | Identifying Y-Chromosome Haplogroups

[![python](
https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue.svg)](
https://docs.python.org)
[![style](
https://img.shields.io/badge/style-black-blue.svg)](
https://github.com/psf/black)
[![imports](
https://img.shields.io/badge/imports-isort-blue.svg)](
https://pycqa.github.io/isort)
[![docs](
https://img.shields.io/badge/docs-pydocstyle-blue.svg)](
https://github.com/PyCQA/pydocstyle)
[![format](
https://img.shields.io/badge/format-flake8-blue.svg)](
https://flake8.pycqa.org/en/latest/index.html)
[![mypy](
https://img.shields.io/badge/annotations-mypy-blue.svg)](
https://github.com/python/mypy)

David Poznik, 23andMe

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
bioRxiv [preprint](http://biorxiv.org/content/early/2016/11/19/088716):

    Poznik GD. 2016. Identifying Y-chromosome haplogroups in arbitrarily large samples
    of sequenced or genotyped men. bioRxiv doi: 10.1101/088716

To learn more about the software, please see the manual, [`yhaplo_manual.pdf`](
./yhaplo_manual.pdf).

For an overiew of command-line options, install the package and run `yhaplo --help`.


## Installation

### Basic installation

To install:
```sh
git clone git@github.com:23andMe/yhaplo.git
cd yhaplo
pip install --editable .
```

To update:
```sh
cd /path/to/yhaplo
git pull  # Update code
pip install --editable .  # Update version number
```

### Optional dependencies

To include optional dependencies for various features:
* `pip install --editable .[vcf]` Enables running on VCF/BCF input
* `pip install --editable .[plot]` Enables tree plotting
* `pip install --editable .[dev]` Includes all optional dependencies,
  as well as development tools (e.g., `pytest`)

To install multiple optional features, use a comma-separated list. For example:
```sh
pip install --editable .[vcf,plot]
```


## Testing

### Running on example data

To run on example text data:
```sh
yhaplo --example_text
```

The `--example_text` option tells `yhaplo` to run on a subset of 1000 Genomes data
in sample-major text format. It also sets the <nobr>`--all_aux_output`</nobr> flag
to produce all auxiliary output.

Similarly, to run on example VCF data:
```sh
yhaplo --example_vcf
```

### Unit tests

To run unit tests:
```sh
make test
```


## Caveats

Please note the following caveats before running `yhaplo`:

* `yhaplo` does not check for sex status; it assumes all individuals are male.
* `yhaplo` expects SNP coordinates consistent with the hg19/GRCh37 reference assembly.
* `yhaplo` expects data at a reasonable number of ISOGG SNPs. This assumption is violated by:
  * Variants-only sequence data
  * Very low-coverage sequencing
  * Genotyping arrays with few Y-chromosome probes

If, for a given individual, `yhaplo` observes no derived alleles at ISOGG SNPs on the upper
branches of the Y-chromosome phylogeny, it will call the individual haplogroup "A,"
since all human Y-chromosome lineages are technically sublineages of A.
Before concluding that the individual sample belongs to paragroup A (which
includes haplogroups A00, A0, A1a, and A1b1), run with the `--anc_snps` option, and check the
auxiliary output for ancestral alleles at haplogroup-BT SNPs. If you do not see any,
your data set probably violates one or more of the assumptions listed above.

In particular, "variants-only" VCF files restrict to SNPs at which alternative alleles
were observed, but ref/alt status is unimportant to `yhaplo`. What is important is
ancestral/derived status. The reference sequence contains many derived alleles,
and `yhaplo` will not be happy if you discard these valuable data. So please emit all
confident sites when calling variants. To limit file size, you could safely restrict to
positions in `output/isogg.snps.unique.DATE.txt`, as these are the only SNPs `yhaplo`
considers. To generate this file, just run `yhaplo` with no arguments.


## Input

The following input file types are supported:
* Indexed BCF: `.bcf`, `.bcf.csi`
* Indexed VCF: `.vcf.gz`, `.vcf.gz.tbi`
* Sample-major text: `.genos.txt` or `.genos.txt.gz`
    * Row 0: Physical coordinates (GRCh37)
    * Column 0: Individual identifiers
    * Cell (*i*, *j*): Genotype for individual *i* at position *j*.<br>
      Values include {"A", "C", "G", "T", "."}, with "." indicating an unobserved value.

In addition, the API supports running on a mapping of individual identifiers to 23andMe ablocks.


## Output

All output file formats are described in detail in [`yhaplo_manual.pdf`](
./yhaplo_manual.pdf).

The two primary output files are:

1. `log.project_name.txt` Log file containing details of the run
2. `haplogroups.project_name.txt` Haplogroup calls. The 4 columns are:
    1. ID
    2. Haplogroup short form, with the name of a SNP observed in the derived state
    3. Haplogroup short form, with the name of a representative SNP
    4. Haplogroup long form, using Y-Chromosome Consortium nomenclature

`yhaplo` also produces a number of SNP tables, tree files, and auxiliary output files.<br>
Please see [`yhaplo_manual.pdf`](./yhaplo_manual.pdf) and `yhaplo --help` for details.


## API

See `yhaplo/api/call_haplogroups.py`.


## CLI

The main command-line entry-point is `yhaplo`.
Additional commands include:
* `yhaplo_convert_to_genos`
* `yhaplo_plot_tree`


## Implementation details

### Package data

#### Tree

The primary structure of the Y-chromosome tree is stored in
`yhaplo/data/tree/y.tree.primary.DATE.nwk`.

#### Variants

Variant metadata are stored in `yhaplo/data/variants/`:
* `isogg.DATE.txt` Phylogenetically informative SNPs scraped directly from the ISOGG website.<br>
  `yhaplo` resolves errors and formatting inconsistencies and emits cleaned versions:
  `isogg.snps.cleaned.DATE.txt`, `isogg.snps.unique.DATE.txt`.<br>
  See [`yhaplo_manual.pdf`](./yhaplo_manual.pdf) for details.
* `isogg.correct.*.txt` Corrections to ISOGG data
* `isogg.multiallelic.txt` Physical coordinates of multiallelic sites to be excluded
* `isogg.omit.*.txt` SNPs to drop due to inconsistencies observed in test data
* `isogg.split.txt` Not currently used
* `preferred.snp_names.txt` List of preferred SNP names
* `representative.SNPs.*.txt` SNPs deemed representative of corresponding haplogroups


### Classes

#### Trees

The `Tree` class is defined in `tree.py`. It:
* Parses a Newick file to build primary tree
* Parses ISOGG table to add SNPs to nodes and grow tree
* Finds the derived path leading from the root to an individual
* Knows root, depth, haplogroup-to-node mappings, etc.

#### Nodes

The `Node` class is defined in `node.py`. It:
* Represents a phylogenetic branch
* Knows parent, children, SNPs, etc.

#### SNPs

The `SNP` class and related classes are defined in `snp.py`:
* `SNP` Knows position, ancestral and derived alleles, node, etc.

#### Samples

The `Sample` class and its subclasses are defined in `sample.py`:
* `Sample` Knows genotypes and haplogroup of an individual
* `TextSample` Subclass for sample-major text input
* `VCFSample` Subclass for VCF/BCF input

#### Paths

The `Path` class is defined in `path.py`. It represents a ath through a tree and stores:
* The next node to visit
* A list of SNPs observed in the derived state
* The most derived SNP observed
* The number of ancestral alleles encountered

#### Configuration

The `Config` class is defined in `config.py`.
It is a container for parameters, command-line options, and filenames.
