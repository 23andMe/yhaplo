# Changelog for `yhaplo`

Format based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)


## [2.1.9]

### Fixed
- `Node.hg_snp` values would drift upon repeated `Tree` instantiation

[2.1.9]: https://github.com/23andMe/yhaplo/compare/2.1.8..2.1.9


## [2.1.8]

### Added
- `--haplogroup_query` option

[2.1.8]: https://github.com/23andMe/yhaplo/compare/2.1.7..2.1.8


## [2.1.7] - 2024-02-23

### Added
- `SNP.info` property, called from `Tree.query_snp_path`

[2.1.7]: https://github.com/23andMe/yhaplo/compare/2.1.6..2.1.7


## [2.1.6] - 2024-02-07

### Added
- Python 3.12 support
- `__str__` and `__repr__` methods
- SNP-based haplogroup in MRCA output

### Changed
- Newick representation is now unrotated by default

### Removed
- `setuptools_scm`-generated version file

### Fixed
- When generating Newick representation, recalculate maximum depth to support pruned trees

[2.1.6]: https://github.com/23andMe/yhaplo/compare/2.1.4..2.1.6


## [2.1.4] - 2024-01-29

### Added
- Support for querying multiple SNPs from the command line
- Support for generating aligned-tip Newick representations of subtrees
- Support for node lookup by SNP-based haplogroup label
- `Node` methods:
  - `Node.iter_depth_first`
  - `Node.iter_breath_first`
  - `Node.remove_children`
- `tree.get_bounded_subtree` function

### Changed
- `Tree` constructor's `Config` parameter is now optional
- `Node.is_root` and `Node.is_leaf` are now properties

### Fixed
- Output format of `--snp_query` option

[2.1.4]: https://github.com/23andMe/yhaplo/compare/2.1.0...2.1.4


## [2.1.0] - 2024-01-12

This release improves haplogroup calling by identifying and correcting various errors in
the ISOGG variant metadata and, internally, by pruning poorly performing v5 SNPs.

### Changed
- Identify and correct name, position, and mutation errors in ISOGG variant metadata

[2.1.0]: https://github.com/23andMe/yhaplo/compare/2.0.2...2.1.0


## [2.0.2] - 2023-09-15

This is a major clean-up and refactoring release.
Core logic has not changed, and output should be equivalent to prior versions.
The key changes from an end-user perspective are BCF support, a cleaner API,
and faster processing of most input types.

### Added
- BCF support
- Automated tests
- Optional dependencies
- `Sample` subclasses: `TextSample`, `VCFSample`
- `CHANGELOG.md`

### Changed
- Lint code
- Update pre-commit hooks
- Update `Makefile` and configuration files
- Refactor for PEP-8 compliance (snake case, etc.)
- Update directory structure
- Modernize packaging and infer version dynamically
- Namespace command-line entry points: `yhaplo`, `yhaplo_convert_to_genos`, `yhaplo_plot_tree`
- Replace static methods
- Clean up logging and use file handlers
- Use f-strings
- Reformat docstrings
- Add type annotations
- Use `importlib.resources` to load metadata files
- Move example input from package to `tests/fixtures/`
- Update `README.md`, `README.23andMe.md`, and `yhaplo_manual.pdf`
- Speed up sample-major file processing
- Use Pysam to process VCF/BCF input

### Removed
- Support for Python 2 and Python 3.8

[2.0.2]: https://github.com/23andMe/yhaplo/compare/1.1.2..2.0.2

