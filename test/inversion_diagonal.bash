#!/usr/bin/env bash

# Model definition
input_data=../INPUT/Recorded_time_sources_checkerboard_lr_100x100_arma.txt
input_matrix=../INPUT/matrix_checkerboard_lr_100x100_arma.txt

# Output file
name=new_alg_mass_diagonal_10e3_checkerboard_100x100
output_samples=../OUTPUT/${name}.txt
output_trajectory=../OUTPUT/${name}_trajecory.txt
output_log=../OUTPUT/${name}.log
output_plot=../OUTPUT/${name}.pdf

# plot (0) or save (1)
pl_or_sv=1

# Prior info
means=0.000667
std_dev=0.001

# Tuning parameters
algorithm_type=1 # 1 for new, 0 for neal
mass_matrix_type=1 # 0 for complete, 1 for diagonal, 2 for unit
ergodicity=1 # 1 for enforcing ergodicity, 0 for not
temperature=450
adapttimestep=1
time_step=nan # nan for default, is overridden by adapttmestep
number_of_samples=$((10**4))


# Run inversion
../bin/hmc_sampler \
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
python ../src/visualization/plotTomography.py ${output_samples} ${pl_or_sv} ${output_plot}
