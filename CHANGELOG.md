# Changelog

## Unreleased

### Added

### Fixed

- Fixed multiple build issues on Mac OS
- Removed invalid unit test

## [1.0.0] - 2023-06-09

### Added

- Added direct QR factorization to solve augmented Jacobian systems
- Extended dynamic functions to include constraints
- Adapted HiGHS to use newly-introduced integer constants
- Refactored sparse matrix
- Added support for LAPACK factorization
- Improved CUTEst reporting
- Added slack bases to initial LP solves
- Restricted derivative checks to first iteration
- Consolidated Params / Options into Settings
- Added notebooks as documentation
- Added capability to read settings from files
- Refactored / simplified python interface
- Removed unused objective dual
- Improved prints in solver

### Fixed

- Fixed error in computation of Newton objective
- Fixed error in preprocessing forcing constraints
- Fixed error in constraint classes in Python interface
- Fixed bug in python "minimize" evaluation
- Fixed bug in python interface settings

## [0.1.0] - 2022-07-07

Initial release
