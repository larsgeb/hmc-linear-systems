#!/usr/bin/env bash

# Model definition
input_data=INPUT/Recorded_time_sources_checkerboard_lr_10x10_arma.txt
input_matrix=INPUT/matrix_checkerboard_lr_10x10_arma.txt

# Output file
output_samples=OUTPUT/samples_10x10.txt
output_plot=PLOTS/samples_10x10.pdf

# Prior info
means=0.000667
std_dev=0.001

# Tuning parameters
algorithm_type=1 # 1 for new, 0 for neal
mass_matrix_type=0 # 0 for complete, 1 for diagonal, 2 for unit
ergodicity=1 # 1 for enforcing ergodicity, 0 for not
temperature=10
number_of_samples=$((10**5))


# Run inversion
./runSampling \
    -im ${input_matrix} \
    -id ${input_data} \
    -os ${output_samples} \
    -ns ${number_of_samples} \
    -t ${temperature} \
    -e ${ergodicity} \
    -means ${means} \
    -std ${std_dev} \
    --massmatrixtype ${mass_matrix_type} \
    -an ${algorithm_type}


# Visualize data
python plotTomography.py ${output_samples}