#include <iostream>
#include <stdio.h>
#include <Stats.hpp>
#include <assert.h>
#include <math.h>
#include <limits>

void precompute_pearson_data(int sample_data_length, const double* __restrict__ sample_data, double* __restrict__ precomputed_pearson_data){

    //Calculate sum and average of data samples
    double sum = 0, avg = 0;
    for(int i = 0; i < sample_data_length; i++)
        sum += sample_data[i];

    avg =  sum / sample_data_length;

    //Calculate standard deviation of data samples
    double factor_sum = 0;
    for(int i = 0; i < sample_data_length; i++)
        factor_sum += pow((sample_data[i] - avg),2);

    double stddev = 0;
    stddev = sqrt(factor_sum/sample_data_length);

    //Precompute pearson data
    for(int i = 0; i < sample_data_length; i++){
        if(fabs(stddev) < 2* std::numeric_limits<double>::min())
            precomputed_pearson_data[i] = 0;
        else
            precomputed_pearson_data[i] = (sample_data[i] - avg)/(stddev * sample_data_length);
    }
}

double pearsoncorr_from_precomputed(int n, const double* __restrict__  v1, const double* __restrict__  v2){
    double sum = 0;
    for(int i = 0; i < n; i++){
        sum += v1[i] * v2[i];
    }
    return sum*n;
}


