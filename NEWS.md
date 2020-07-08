# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.0] - 2020-07-08

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
