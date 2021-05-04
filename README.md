# SLEQP

## Dependencies

For the solution of the linear programming problem, one of the following packages is required:

* [Gurobi](https://www.gurobi.com/)
* [SoPlex](https://soplex.zib.de/)

For factorizations, one of the following libraries is required:

* [UMFPACK](http://faculty.cse.tamu.edu/davis/suitesparse.html)
* [MUMPS](http://mumps.enseeiht.fr/)
* [MA27](https://www.hsl.rl.ac.uk/archive/specs/ma27.pdf)
* [MA57](https://www.hsl.rl.ac.uk/archive/specs/ma57.pdf)
* [MA86](https://www.hsl.rl.ac.uk/ipopt/)
* [MA97](https://www.hsl.rl.ac.uk/ipopt/)

Additional dependencies:
* [trlib](https://github.com/felixlen/trlib) for nonconvex trust region QP solving

Optional dependencies for the python bindings:
* [Python](https://www.python.org/) itself, including `setuptools`, `numpy`, and `scipy`
* [Cython](https://cython.org/) to build the extension

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
* `SLEQP_LPS`: Set to specify a linear programming solver
* `SLEQP_FACT`: Set to specify a factorization library

## References

* Byrd, Richard H., et al. "An algorithm for nonlinear optimization using linear programming and equality constrained subproblems." Mathematical Programming 100.1 (2003): 27-48.
* Lenders, Felix, Christian Kirches, and Andreas Potschka. "trlib: A vector-free implementation of the GLTR method for iterative solution of the trust region problem." Optimization Methods and Software 33.3 (2018): 420-449.
