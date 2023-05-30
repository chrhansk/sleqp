# Changelog

## [Unreleased]

### Added

- Added direct QR factorization to solve augmented Jacobian systems
- Extended dynamic functions to include constraints
- Adapted HiGHS to use newly-introduced integer constants
- Refactored sparse matrix
- Added support for LAPACK factorization
- Improved CUTEst reporting
- Added slack bases to initial LP solves

### Fixed

- Fixed error in computation of Newton objective
- Fixed error in preprocessing forcing constraints
- Fixed error in constraint classes in Python interface

## [0.1.0] - 2022-07-07

Initial release
