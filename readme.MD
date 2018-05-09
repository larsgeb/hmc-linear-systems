
## Hamiltonian Monte Carlo - Probabilistic inversion of linear models

This program is built to incorporate MCMC inversion using Hamiltonian Mechanics. Essentially, any amount of
parameters and any linear forward model matrix can be supplied, as long as one can describe the measurement errors, data
uncertainty and misfit function.

The latest implementation (on the 'armadillo' branch) uses the Armadillo Linear Algebra Library 
to do multithread linear algebra operations. This code has been succesfully tested for models
up to 10201 parameters. A detailed documentation is as of yet not complete, but for a general idea how
to set up an inversion i do recommend looking in test.

### Current work
The effectiveness of the method depends mostly on the choice of mass matrix. It turns out that the
optimal mass matrix is the systems Hessian. For a linear problem, solving the Hessian of the Bayesian 
formulation is the same as obtaining the full solution. It is therefore quite backwards to solve the
linear system by computing the Hessian, after which we can use this as a mass matrix to compute the 
solution.

An alternative is having an estimate of the diagonal entries of the Hessian; the variances of each 
model parameter. These are relatively cheap to compute if the prior and data covariances are diagonal. 
This is currently being implemented.

### Prerequisites

Make sure that BLAS and LaPACK are available on your system, as well as the Armadillo library. 

### Compilation

All build targets can easily be prepared with CMake. Using at least CMake 3.7 allows you to run (from the
root directory);

```{r, engine='bash', Running CMake}
cmake .
```
This prepares the platform dependent compilation instructions. For example, on Linux one would now run for
 all separate targets:

 ```{r, engine='bash', Building targets}
 make <target>
 ```

The targets can be found in CMakeLists.txt, but are respectively the sampler (hmc_sampler) and generate 
synthetics (gen_syn).

### Execution

All targets are compiled into bin/, and should be run from there too. For ease of use, several bash scripts
 are included in src/test/, which should again be run from bin. All visualization scripts should 
also be run from this folder.