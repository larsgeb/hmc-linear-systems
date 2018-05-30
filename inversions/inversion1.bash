#!/usr/bin/env bash

# Model definition
input_A=A.txt
input_B=B.txt
input_C=C.txt

# Output file
name=inversion1
output_samples=${name}/samples.txt
output_trajectory=${name}/trajectory.txt
output_log=${name}/${name}.log

# Tuning parameters
mass_matrix_type=0 # 0 for complete, 1 for diagonal, 2 for unit
temperature=1
adapt_time_step=0
time_step=0.1 # nan for default, is overridden by adapttmestep
number_of_samples=$((1000))

# Run inversion
./hmc_sampler \
    -nt 50 \
    -ia ${input_A} \
    -ib ${input_B} \
    -ic ${input_C} \
    -os ${output_samples} \
    -ot ${output_trajectory} \
    -ns ${number_of_samples} \
    -t ${temperature} \
    -dt ${time_step} \
    -means ${means} \
    -at ${adapt_time_step} \
    -std ${std_dev} \
    --massmatrixtype ${mass_matrix_type} \
    2>&1 | tee ${output_log}

sed -i 's/\x1b\[[0-9;]*m//g' ${output_log}
sed -i '/\[/d' ${output_log}