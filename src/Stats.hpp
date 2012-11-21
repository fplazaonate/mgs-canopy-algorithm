#ifndef STATS
#define STATS


double pearsoncorr_from_precomputed(int n, const double* v1, const double* v2);

double morten_pearsoncorr(int n, const double* v1, const double* v2);

double* precompute_pearson_data(const double* sample_data);
#endif
