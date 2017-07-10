
### HMC Inversion of simple (synthetic) VSP data

This set of programs are built to incorporate MCMC inversion using Hamiltonian Mechanics. Essentially, any amount of
parameters and any forward model can be supplied, as long as one can descirbe the measurement errors, data uncertainty
and misfit function. This could be used as a framework for building any inversion using HMC.

## Install

First generate data using the make_data.cpp file. Compile it against the aux class and you're good to go.

The inversion is done using the sampling.cpp file, which should be compiled against the aux and montecarlo class. These are not
finished yet.

## Settings

The geometry setup of the receivers as well as the initial starting model/actual model used for synthetics can all be
edited in the INPUT folder. As of yet, the synthetics model is the same as the starting model.
