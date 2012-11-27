#ifndef STATS_TEST
#define STATS_TEST

#include <stdio.h>
#include <time.h>

#include <boost/test/unit_test.hpp>

#include <Stats.hpp>


using namespace boost::unit_test;

void fill_array_with_random_data(double* array, int len){
    for(int i =0; i < len; i++){
        array[i] = ((double)(rand() % 10000))/100;
    }
}

void copy_array_values(double* source, double* dest, int len){
    for(int i =0; i < len; i++){
        dest[i] = source[i];
    }
}

void print_array(double* arr, int len){
    std::cout << "<<<<<" << std::endl;
    for(int i =0; i < len; i++){
        std::cout << arr[i] << "\t";
        if(!len%10)
            std::cout << std::endl;

    }
    std::cout << "<<<<<" << std::endl;
}

BOOST_AUTO_TEST_CASE( test_pearson_precomputation_simple ){
    int num_samples = 100;
    srand( time(NULL) );

    double* data1 = new double[num_samples];
    double* data2 = new double[num_samples];
    double* pc_data1;
    double* pc_data2;
    double* pc_data_morten_1;
    double* pc_data_morten_2;

    double morten_pearson;
    double precomp_pearson;
    double precomp_morten_pearson;
    double avg_precomp_pearson = 0;

    //
    //Perfect correlation
    //
    for(int i = 0; i < num_samples; i++){
        data1[i] = i;
        data2[i] = 2*i;
    }

    pc_data1 = precompute_pearson_data(num_samples, data1);
    pc_data2 = precompute_pearson_data(num_samples, data2);

    pc_data_morten_1 = precompute_pearson_data_morten(num_samples, data1);
    pc_data_morten_2 = precompute_pearson_data_morten(num_samples, data2);

    precomp_pearson = pearsoncorr_from_precomputed(num_samples, pc_data1, pc_data2);
    precomp_morten_pearson = pearsoncorr_from_precomputed_morten(num_samples, pc_data_morten_1, pc_data_morten_2);
    morten_pearson = morten_pearsoncorr(num_samples, data1, data2);

    BOOST_CHECK_CLOSE( precomp_pearson, 1, 1);
    BOOST_CHECK_CLOSE( precomp_morten_pearson, 1, 1 );
    BOOST_CHECK_CLOSE( morten_pearson, 1, 1 );

    delete pc_data1;
    delete pc_data2;

    delete pc_data_morten_1;
    delete pc_data_morten_2;

    //
    //Perfect negative correlation
    //
    for(int i = 0; i < num_samples; i++){
        data1[i] = i;
        data2[i] = (-2)*i;
    }

    pc_data1 = precompute_pearson_data(num_samples, data1);
    pc_data2 = precompute_pearson_data(num_samples, data2);

    pc_data_morten_1 = precompute_pearson_data_morten(num_samples, data1);
    pc_data_morten_2 = precompute_pearson_data_morten(num_samples, data2);

    precomp_pearson = pearsoncorr_from_precomputed(num_samples, pc_data1, pc_data2);
    precomp_morten_pearson = pearsoncorr_from_precomputed_morten(num_samples, pc_data_morten_1, pc_data_morten_2);
    morten_pearson = morten_pearsoncorr(num_samples, data1, data2);

    delete pc_data1;
    delete pc_data2;

    delete pc_data_morten_1;
    delete pc_data_morten_2;

    BOOST_CHECK_CLOSE( precomp_pearson, -1, 1);
    BOOST_CHECK_CLOSE( precomp_morten_pearson, -1, 1 );
    BOOST_CHECK_CLOSE( morten_pearson, -1, 1 );

    delete data1;
    delete data2;

}

BOOST_AUTO_TEST_CASE( test_pearson_precomputation_simple_random ){

    int num_iterations = 1000; 
    int num_samples = 1000;
    srand( time(NULL) );

    double* data1 = new double[num_samples];
    double* data2 = new double[num_samples];

    double precomp_pearson;
    double avg_precomp_pearson = 0;

    for(int i = 0; i < num_iterations; i++){

        fill_array_with_random_data(data1, num_samples);
        fill_array_with_random_data(data2, num_samples);

        double* pc_data1 = precompute_pearson_data(num_samples, data1);
        double* pc_data2 = precompute_pearson_data(num_samples, data2);

        precomp_pearson = pearsoncorr_from_precomputed(num_samples, pc_data1, pc_data2);

        avg_precomp_pearson+=precomp_pearson;


        delete pc_data1;
        delete pc_data2;
    }
    avg_precomp_pearson/=num_iterations;

    delete data1;
    delete data2;

    BOOST_CHECK_SMALL(avg_precomp_pearson, 0.001);

}

BOOST_AUTO_TEST_CASE( test_pearson_precomputation_against_mortens )
{

    int num_iterations = 1000; 
    int num_samples = 100;
    srand( time(NULL) );

    double* data1 = new double[num_samples];
    double* data2 = new double[num_samples];

    double precomp_pearson;
    double morten_pearson;
    double precomp_morten_pearson;

    double avg_precomp_pearson = 0;
    double avg_morten_pearson = 0;
    
    for(int i = 0; i < num_iterations; i++){

        fill_array_with_random_data(data1, num_samples);
        fill_array_with_random_data(data2, num_samples);

        double* pc_data1 = precompute_pearson_data(num_samples, data1);
        double* pc_data2 = precompute_pearson_data(num_samples, data2);

        double* pc_data_morten_1 = precompute_pearson_data_morten(num_samples, data1);
        double* pc_data_morten_2 = precompute_pearson_data_morten(num_samples, data2);
        
        precomp_pearson = pearsoncorr_from_precomputed(num_samples, pc_data1, pc_data2);
        morten_pearson = morten_pearsoncorr(num_samples, data1, data2);
        precomp_morten_pearson = pearsoncorr_from_precomputed_morten(num_samples, pc_data_morten_1, pc_data_morten_2);

        avg_precomp_pearson+=precomp_pearson;
        avg_morten_pearson+=morten_pearson;

        BOOST_CHECK_CLOSE( precomp_morten_pearson, morten_pearson, 1 );
        BOOST_CHECK_CLOSE( precomp_pearson, morten_pearson, 1 );


        delete pc_data1;
        delete pc_data2;

        delete pc_data_morten_1;
        delete pc_data_morten_2;
    }

    BOOST_CHECK_CLOSE( avg_precomp_pearson, avg_morten_pearson, 1 );

    delete data1;
    delete data2;
}



#endif
