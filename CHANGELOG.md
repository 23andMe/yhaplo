# Changelog for `yhaplo`

Format based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)


## Planned

- Correct ISOGG polarization errors for a few dozen SNPs.


## [Unreleased]

No unreleased changes

[Unreleased]: https://github.com/23andMe/yhaplo/compare/2.0.2...HEAD


## [2.0.2] - 2023-09-15

This is a major clean-up and refactoring release.
Core logic has not changed, and output should be equivalent to prior versions.
The key changes from an end-user perspective are BCF support, a cleaner API,
and faster processing of most input types.

### Added
- BCF support
- Automated tests
- Optional dependencies
- `Sample` subclasses: `TextSample`, `VCFSample`, `AblockSample`
- API for processing ablocks (23andMe internal)
- `Dockerfile` defining image for Batch computes (23andMe internal)
- Compute flow (23andMe internal)
- Script for copying and altering files for open sourcing (23andMe internal)
- `CHANGELOG.md`

### Changed
- Lint and update pre-commit hooks
- Set up Drone CI (23andMe internal)
- Set up `tox` testing (23andMe internal)
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
- Speed up ablock processing (23andMe internal)
- Use Pysam to process VCF/BCF input
- Map physical coordinates to block indexes (23andMe internal)
- Handle platform SNPs natively (23andMe internal)

### Removed
- Support for Python 2 and Python 3.8
- Use of research-environment utilities (23andMe internal)

[2.0.2]: https://github.com/23andMe/yhaplo/compare/1.1.2..2.0.2

