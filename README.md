# SLEQP

## Introduction

SLEQP is a software package for
large-scale [nonlinear optimization](https://en.wikipedia.org/wiki/Nonlinear_programming).
It is designed to find (local) solutions of mathematical optimization problems of the form

```
   min     f(x)
  x ∈ ℝⁿ

s.t.         l ≤ c(x) ≤ u
           x_l ≤  x   ≤ x_u
```

where ```f: ℝⁿ --> ℝ``` is the optimization objective, and ```c: ℝⁿ --> ℝᵐ```
are optimization constraints.

The vectors `l, u ∈ ℝᵐ` denote the lower and upper bounds on the
variables, while the vectors `x_l, x_u ∈ ℝⁿ` are bounds on the
nonlinear constraints. The functions `f(x)` and `c(x)` can be
nonlinear and nonconvex, but should be twice continuously
differentiable.

## Dependencies

For the solution of the linear programming problem, one of the
following linear programming libraries is required:

* [Gurobi](https://www.gurobi.com/)
* [HiGHS](https://www.highs.dev/)
* [SoPlex](https://soplex.zib.de/)

For factorizations, one of the following factorization libraries is required:

* [Umfpack](http://faculty.cse.tamu.edu/davis/suitesparse.html)
* [MUMPS](http://mumps.enseeiht.fr/)
* [MA27](https://www.hsl.rl.ac.uk/archive/specs/ma27.pdf)
* [MA57](https://www.hsl.rl.ac.uk/archive/specs/ma57.pdf)
* [MA86](https://www.hsl.rl.ac.uk/ipopt/)
* [MA97](https://www.hsl.rl.ac.uk/ipopt/)
* [SPQR](http://faculty.cse.tamu.edu/davis/suitesparse.html)

Additional dependencies:
* [trlib](https://github.com/felixlen/trlib) for nonconvex trust region QP solving

Optional dependencies for the python bindings:
* [Python](https://www.python.org/) itself, including `setuptools`, `numpy`, and `scipy`
* [Cython](https://cython.org/) to build the extension

Optional dependencies for the mex bindings:
* [MATLAB](https://www.mathworks.com/products/matlab.html), *or*
* [Octave](https://www.gnu.org/software/octave/index)

Optional dependencies for the AMPL bindings:
* The AMPL solver library [ASL](https://github.com/coin-or-tools/ThirdParty-ASL)

Optional dependencies for the unit tests:
*  [Check](https://libcheck.github.io/check/)

## Installation

In order to compile this package, use the following sequence of commands:

```
mkdir build && cd build
cmake .. <OPTIONS>
make
[make build_tests && make test]
```

### Options

Use the following options to customize the build process:

* `SLEQP_ENABLE_UNIT_TESTS`: Enables the unit tests (default: `On`)
* `SLEQP_ENABLE_CUTEST`: Enables the CUTest suite (default: `Off`)
* `SLEQP_ENABLE_PYTHON`: Enables the build of the python bindings (default : `On`)
* `SLEQP_LPS`: Set to specify a linear programming library
* `SLEQP_FACT`: Set to specify a factorization library
* `SLEQP_ENABLE_MATLAB_MEX`: Enables the build of mex bindings using MATLAB (default : `Off`)
* `SLEQP_ENABLE_OCTAVE_MEX`: Enables the build of mex bindings using Octave (default : `Off`)
* `SLEQP_ENABLE_AMPL`: Enables the build of the AMPL interface (default: `Off`)

## References

* Byrd, Richard H., et al. "An algorithm for nonlinear optimization using linear programming and equality constrained subproblems." Mathematical Programming 100.1 (2003): 27-48.
* Lenders, Felix, Christian Kirches, and Andreas Potschka. "trlib: A vector-free implementation of the GLTR method for iterative solution of the trust region problem." Optimization Methods and Software 33.3 (2018): 420-449.
