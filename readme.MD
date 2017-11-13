
## hmc Inversion of simple linear models

This set of programs is built to incorporate MCMC inversion using Hamiltonian Mechanics. Essentially, any amount of
parameters and any linear forward model matrix can be supplied, as long as one can describe the measurement errors, data
uncertainty and misfit function.

The latest implementation (on the 'armadillo' branch) uses the Armadillo Linear Algebra Library 
to do multithread linear algebra operations. This code has been succesfully tested for models
up to 10201 parameters. A detailed documentation is as of yet not complete (see docs/), but looking 
through runSampling.cpp one will encounter most settings as member functions of the settings object.

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

There are two build targets:
1. generateSynthetics (deprecated).

   This program will read the forward matrix in INPUT/matrix.txt and generate synthetics which are saved in 
   synthetics.txt. Model parameters are as of yet hardcoded. NOTE: This code must be adapted still for Armadillo.

2. runSampling.

   This program will take samples from the posterior distribution described by hardcoded 
   prior information, the forward model from a text file and synthetic data 
   from another text files. (with hardcoded data-covariance). The output are the samples in 
   OUTPUT/samples.txt and the trajectory of the last sample in OUTPUT/trajectory.txt. 
   Some settings and diagnostics are written to the console.

### Code documentation

Documentation generated using doxygen is found in the docs/ subdirectory, but is also found hosted at [my GitHub 
Pages](www.larsgeb.github.io).
