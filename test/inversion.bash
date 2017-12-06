#!/usr/bin/env bash

# Model definition
input_data=../INPUT/Recorded_time_sources_checkerboard_lr_10x10_arma.txt
input_matrix=../INPUT/matrix_checkerboard_lr_10x10_arma.txt

# Output file
name=checkerboard_10x10
output_samples=../OUTPUT/${name}.txt
output_trajectory=../OUTPUT/${name}_trajecory.txt
output_log=../OUTPUT/${name}.log
output_plot=../OUTPUT/${name}.pdf

# Prior info
means=0.000667
std_dev=0.001

# Tuning parameters
algorithm_type=1 # 1 for new, 0 for neal
mass_matrix_type=0 # 0 for complete, 1 for diagonal, 2 for unit
ergodicity=1 # 1 for enforcing ergodicity, 0 for not
temperature=10
number_of_samples=$((10**4))


# Run inversion
../bin/runSampling \
    -im ${input_matrix} \
    -id ${input_data} \
    -os ${output_samples} \
    -ot ${output_trajectory} \
    -ns ${number_of_samples} \
    -t ${temperature} \
    -e ${ergodicity} \
    -means ${means} \
    -std ${std_dev} \
    --massmatrixtype ${mass_matrix_type} \
    -an ${algorithm_type} \
    2>&1 | tee ${output_log}

# Visualize data
python ../src/visualization/plotTomography.py ${output_samples}
