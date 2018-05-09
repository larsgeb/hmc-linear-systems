#!/usr/bin/env bash

# Model definition
input_data=observed_data.bin
input_matrix=model_matrix.bin

# Output file
name=inversion_1
output_samples=OUTPUT/${name}.txt
output_trajectory=OUTPUT/${name}_trajecory.txt
output_log=OUTPUT/${name}.log

# Prior info
means=1
std_dev=10

# Tuning parameters
algorithm_type=0 # 1 for new, 0 for neal
mass_matrix_type=1 # 0 for complete, 1 for diagonal, 2 for unit
ergodicity=1 # 1 for enforcing ergodicity, 0 for not
temperature=1
adapttimestep=0
time_step=0.5 # nan for default, is overridden by adapttmestep
number_of_samples=$((10**3))

# Run inversion
../bin/hmc_sampler \
    -nt 100 \
    -im ${input_matrix} \
    -id ${input_data} \
    -os ${output_samples} \
    -ot ${output_trajectory} \
    -ns ${number_of_samples} \
    -t ${temperature} \
    -e ${ergodicity} \
    -dt ${time_step} \
    -means ${means} \
    -at ${adapttimestep} \
    -std ${std_dev} \
    --massmatrixtype ${mass_matrix_type} \
    -an ${algorithm_type} \
    2>&1 | tee ${output_log}

# Visualize data
python ../src/visualization/plotPosterior.py ${output_samples}