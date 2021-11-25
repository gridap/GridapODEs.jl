# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Changed

- Hiding the creation of `TransientFESolver` from user code. Since PR [#59](https://github.com/gridap/GridapODEs.jl/pull/59).

## [0.7.0] - 2021-09-21

### Added
  - Added Newmark method in PR [#55](https://github.com/gridap/GridapODEs.jl/pull/55)
  - Updated to Gridap v0.16 in PR [#54](https://github.com/gridap/GridapODEs.jl/pull/54)

## [0.6.0] - 2020-12-15

### Changed
  - Updated to Gridap v0.15 in PR [#43](https://github.com/gridap/GridapODEs.jl/pull/43)
  
## [0.5.0] - 2020-08-24

### Fixed

  - Append `matdata` from `jacobian` and `jacobian_t` to get the right sparsity pattern. Since PR [#21](https://github.com/gridap/GridapODEs.jl/pull/21)
  - Bug-fixed skeleton terms. Since PR [#24](https://github.com/gridap/GridapODEs.jl/pull/24)

### Changed

  - The `MultiFieldFESpace` method that takes a vector of spaces that include at least one `TransientFESpace` is called `TransientMultiFieldFESpace`. Since PR [#27](https://github.com/gridap/GridapODEs.jl/pull/27)
  - Updated to Gridap v0.13.

## [0.4.0] - 2020-07-08

### Added

  - Support for computation of partial time derivatives since PR [#20](https://github.com/gridap/GridapODEs.jl/pull/20)
  - Updated to Gridap v0.12.0

## [0.3.0] - 2020-06-29

### Added

  - New constant and affine version of operators via a new Trait `OperatorType` being used in `ODEOperator`, `ODEOpFromFEOp` and `TransientFEOperator`. Since PR [#18](https://github.com/gridap/GridapODEs.jl/pull/18)

## [0.2.1] - 2020-06-16

### Changed

  - Tests updated for new Gridap inner product format

## [0.2.0] - 2020-04-26

### Added

  - Support to automatic differentiation
  - Time integration of multivalue problems
  - Implementation of a ForwardEuler solver

### Changed

### Fixed

  - Bug in ThetaMethod

## [0.1.0] - 2020-04-18

It is the first fully functional version.
