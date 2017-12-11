
## Hamiltonian Monte Carlo - Probabilistic inversion of linear models

This program is built to incorporate MCMC inversion using Hamiltonian Mechanics. Essentially, any amount of
parameters and any linear forward model matrix can be supplied, as long as one can describe the measurement errors, data
uncertainty and misfit function.

The latest implementation (on the 'armadillo' branch) uses the Armadillo Linear Algebra Library 
to do multithread linear algebra operations. This code has been succesfully tested for models
up to 10201 parameters. A detailed documentation is as of yet not complete, but for a general idea how
to set up an inversion i do recommend looking in test/.

### Prerequisites

Make sure that BLAS and LaPACK are available on your system, as well as the Armadillo library. Look
at the branch 'euler' how to use Armadillo of one has the libraries only in the project.

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

There is one target as of yet:

1. runSampling

   This program will take samples from the posterior distribution described by supplied 
   prior information, the forward model from a text file and synthetic data 
   from another text file (with hardcoded data-covariance). The output are the samples in 
   OUTPUT/samples.txt and the trajectory of the last sample in OUTPUT/trajectory.txt. 
   Some settings and diagnostics are written to the console.

### Code documentation

Documentation generated using doxygen is found in the docs/ subdirectory, but is also found hosted at [my GitHub 
Pages](www.larsgeb.github.io).
