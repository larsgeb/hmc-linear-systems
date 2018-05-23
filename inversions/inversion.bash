#!/usr/bin/env bash

# Model definition
input_A=A
input_B=B
input_C=C

# Output file
name=inversion_1
output_samples=${name}.txt
output_trajectory=${name}_trajectory.txt
output_log=${name}.log

# Prior info
means=1
std_dev=10

# Tuning parameters
algorithm_type=0 # 1 for new, 0 for neal
mass_matrix_type=1 # 0 for complete, 1 for diagonal, 2 for unit
temperature=1
adapttimestep=1
time_step=nan # nan for default, is overridden by adapttmestep
number_of_samples=$((1000))

# Run inversion
../bin/hmc_sampler \
    -nt 10 \
    -ia ${input_A} \
    -ib ${input_B} \
    -ic ${input_C} \
    -os ${output_samples} \
    -ot ${output_trajectory} \
    -ns ${number_of_samples} \
    -t ${temperature} \
    -dt ${time_step} \
    -means ${means} \
    -at ${adapttimestep} \
    -std ${std_dev} \
    --massmatrixtype ${mass_matrix_type} \
    -an ${algorithm_type} \
    2>&1 | tee ${output_log}

