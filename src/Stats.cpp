#include <Stats.hpp>
#include <assert.h>
#include <math.h>

double* precompute_pearson_data(const double* sample_data){

    //Calculate & verify sample_data length, Allocate new array of the same size
    int sample_data_length = (int)(sizeof(sample_data)/sizeof(double));
    //assert(sample_data_length >= 2);
    double* precomputed_pearson_data = new double[sample_data_length];

    //Calculate sum and average of data samples
    double sum = 0, avg = 0;
    for(int i = 0; i < sample_data_length; i++)
        sum += sample_data[i];

    avg /= sample_data_length;

    //Calculate standard deviation of data samples

    double factor_sum = 0;
    for(int i = 0; i < sample_data_length; i++)
        factor_sum += pow((sample_data[i] - avg),2);

    double stddev = 0;
    stddev = sqrt(factor_sum/sample_data_length);

    //Precompute pearson data
    for(int i = 0; i < sample_data_length; i++)
        precomputed_pearson_data[i] = (sample_data[i] - avg)/(stddev * (sample_data_length - 1));

    return precomputed_pearson_data;
}

double pearsoncorr_from_precomputed(int n, const double* v1, const double* v2){
    double sum = 0;
    for(int i = 0; i < n; i++){
        sum += v1[i] * v2[i];
    }
    return sum;
}



double morten_pearsoncorr(int n, const double* v1, const double* v2){
	int	i;
	double avg_v1=0,avg_v2=0;
	double t, nx, ny;
	double c;

    for(int i = 0; i < n; i++){
        avg_v1+=v1[i];
        avg_v2+=v2[i];
    }
    avg_v1/=n;
    avg_v2/=n;

	t = nx = ny = 0.0;

	for ( i=0;i<n;i++ ) {
		t += ( v1[i] - avg_v1 ) * ( v2[i] - avg_v2 );
		nx += ( v1[i] - avg_v1 ) * ( v1[i] - avg_v1 );
		ny += ( v2[i] - avg_v2 ) * ( v2[i] - avg_v2 );
	}

	if ( nx * ny == 0.0 )
		c = 0.0;
	else
		c = t/sqrt(nx*ny);

	return( c );
}



