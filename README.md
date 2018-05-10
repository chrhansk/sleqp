# SLEQP

### Installation of Packaged Dependencies

`sudo apt install build-essential cmake ccmake libadolc-dev libcolpack-dev libsuitesparse-dev liblapack-dev libblas-dev gfortran python-dev python-numpy python-scipy cython`

Many of the external dependencies need to be downloaded and installed manually.


### Installation of Non-Packaged External Dependencies

For linear program solving, one or more of the following packages are required:
* CLP
* CPLEX
* GuRoBi
* SoPLEX

For nonconvex trust region QP solving, one or more of the following packages are required:
* trlib
* GLTR
If none are available, only trust region QPs will be convexified.

For linear system solving, one or more of the following packages are required:
* LAPACK (dense) (packaged)
* MA27
* MA57
* MA97 
* UMFPACK (packaged)
* MeTiS (optional)

See below for details on the installation procedures of individual packages.


### Installation

```
git clone git@github.com/felixlen/SLEQP.git
mkdir build && cd build
cmake ..
ccmake .
make -j<NUM_OF_CORES>
make install
```

### Development Hints

`git` makes it very easy to have branches and check against master. Goal is always to have a working `master` branch. If you do something, create a new branch by

`git checkout -b mybranch`

and work on that.
You can frequently commit and push to `github`.
Once you're done, open a pull request on `github` to merge your branch into `master`.
Creating the pull request automatically tests it and shows conflict.
If you have the confirmation that everything is fine, merge the pull request in github.
You can automatically delete the create branch after the pull request is merged.

### References

* 
*
*
*

### Non-Packaged External Dependencies

#### LP/GuRoBi
You need to get GuRoBi at gurobi.com. Set the variable `GUROBIDIR` in the pysleqp
`Makefile` to the GuRoBi folder that containts include and lib as subfolders,
for example

    GUROBIDIR = /opt/gurobi650/linux

and make shure that `$GUROBIDIR/lib` is in the `LD_LIBRARY_PATH` environment variable
so that `libgurobi*.so` can be found.

#### LP/CLP

todo.

#### LP/SoPLEX

todo.

#### LP/CPLEX

todo.

#### LS/MeTiS
METIS is a requirement of MA57, which only works with METIS versions < 5.
pysleqp Makefile will automatically download and build METIS for you.
Please acknowledge METIS and its licence at [Karypis Lab](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview).

#### LS/MA27

#### LS/MA57
Get HSL MA57 from the [HSL MA57 webpage](http://www.hsl.rl.ac.uk/catalogue/hsl_ma57.html)
Make sure to obtain HSL MA57 (Fortran 90) and not MA57 (Fortran 77).
Either unpack the HSL MA57 archive to `ThirdParty/hsl_ma57` or copy
the following files from the HSL MA57 archive to `ThirdParty/hsl_ma57/src`:

* `ddeps90.f90`
* `ddeps.f`
* `hsl_ma57d.f90`

#### LS/MA97

todo.

#### TR/trlib

todo.

#### TR/GLTR

Get GALAHAD at the [GALAHAD webpage](http://www.galahad.rl.ac.uk).
Either unpack the GALAHAD archive to `ThirdParty/galahad` or copy
the following files from the GALAHAD archive to `ThirdParty/galahad`,
make sure to preserve relative paths:

* `seds/double.sed`
* `seds/ma57v4.sed`
* `src/auxiliary/norms.f90`
* `src/gltr/gltr.f90`
* `src/lapack/blas_interface.f90`
* `src/lapack/lapack_interface.f90`
* `src/rand/rand.f90`
* `src/roots/roots.f90`
* `src/sort/sort.f90`
* `src/space/space.f90`
* `src/spec/specfile.f90`
* `src/sym/symbols.f90`
