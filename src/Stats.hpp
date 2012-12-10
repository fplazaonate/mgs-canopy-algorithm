#ifndef STATS
#define STATS


double pearsoncorr_from_precomputed(int n, const double* v1, const double* v2);

double morten_pearsoncorr(int n, const double* v1, const double* v2);

void precompute_pearson_data(int sample_data_length, const double* sample_data, double* precomputed_pearson_data);

double* precompute_pearson_data_morten(int sample_data_length, const double* sample_data);

double pearsoncorr_from_precomputed_morten(int n, const double* v1, const double* v2);

#endif
