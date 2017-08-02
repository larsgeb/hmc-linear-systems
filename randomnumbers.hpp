//
// Created by Lars Gebraad on 7/10/17.
//

#ifndef HMC_VSP_RANDOMNUMBERS_HPP
#define HMC_VSP_RANDOMNUMBERS_HPP

const double PI = 3.14159265358979323846264338327;

void randn(double mean, double stdv, double *x1, double *x2);

double randn(double mean, double stdv);

double randf(double min, double max);

#endif //HMC_VSP_RANDOMNUMBERS_HPP
