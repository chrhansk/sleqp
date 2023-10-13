# Changelog

## Unreleased

### Added

- Added function to set enum value from strings in settings
- Added support for coinmumps, absense of MPI
- Added nodiscard attribute to return codes
- Added reference to published PAMM paper

### Fixed

- Added missing symbol export / installation to Windows installation
- Fixed compile / test errors in python extension
- Improved const-correctness of interface functions

## [1.0.1] - 2023-06-26

### Added

- Added feature summary
- Added support for usage via FetchContent()
- Added more documentation of python interface
- Made numpy arrays passed to function callbacks read-only
- Improved restoration phase

### Fixed

- Fixed multiple build issues on Mac OS
- Removed invalid unit test
- Removed source directory variables from cmake-related files
- Fixed error in restoration mode

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
