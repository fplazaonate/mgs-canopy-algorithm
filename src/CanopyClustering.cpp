#include <iostream>
#include <fstream>
#include <ctime>

#include <algorithm>

#include <omp.h>

#include <boost/foreach.hpp>


#include <CanopyClustering.hpp>
#include <Log.hpp>

#include <prog_bar_misc.hpp>

Canopy* CanopyClusteringAlg::create_canopy(Point* origin, vector<Point*>& points, vector<Point*>& close_points, double max_neighbour_dist, double max_close_dist, bool set_close_points){

    std::vector<Point*> neighbours;

    if(set_close_points){
        //Go through all points and set the close points to contain the ones that are "close"
        close_points.clear();//Will not reallocate
        close_points.push_back(origin);
        for(int i=0; i<points.size(); i++){

            Point* potential_neighbour = points[i];
            double dist = get_distance_between_points(origin, potential_neighbour);

            if(dist < max_close_dist && origin != potential_neighbour){

                close_points.push_back(potential_neighbour);

                if(dist < max_neighbour_dist){

                    neighbours.push_back(potential_neighbour);

                }
            } 
        }
    } else {

        for(int i=0; i<close_points.size(); i++){

            Point* potential_neighbour = close_points[i];
            double dist = get_distance_between_points(origin, potential_neighbour);

            //TODO: get rid of this hack
            if((dist < max_neighbour_dist) && (!potential_neighbour->id.compare("!GENERATED!"))){
                neighbours.push_back(potential_neighbour);
            }
        }

    }

    neighbours.push_back(origin);

    Point* center;

    if(neighbours.size() == 1)
        center = origin;
    else
        center = get_centroid_of_points(neighbours);

    return new Canopy(origin, center, neighbours);

}

void CanopyClusteringAlg::filter_clusters_by_single_point_skew(double max_single_data_point_proportion, std::vector<Canopy*>& canopies_to_filter){

    vector<int> canopy_indexes_to_remove;

    for(int i=0; i < canopies_to_filter.size(); i++){
        Point* ccenter = canopies_to_filter[i]->center;
        if(! ccenter->check_if_single_point_proportion_is_smaller_than(max_single_data_point_proportion) )
            canopy_indexes_to_remove.push_back(i);
    }

    std::sort(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());
    std::reverse(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());

    for(int i=0; i < canopy_indexes_to_remove.size(); i++)
        canopies_to_filter.erase(canopies_to_filter.begin() + canopy_indexes_to_remove[i]);

}

void CanopyClusteringAlg::filter_clusters_by_zero_medians(int min_num_non_zero_medians, std::vector<Canopy*>& canopies_to_filter){

    vector<int> canopy_indexes_to_remove;

    for(int i=0; i < canopies_to_filter.size(); i++){
        Point* ccenter = canopies_to_filter[i]->center;
        if(! ccenter->check_if_num_non_zero_samples_is_greater_than_x(min_num_non_zero_medians) )
            canopy_indexes_to_remove.push_back(i);
    }

    std::sort(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());
    std::reverse(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());

    for(int i=0; i < canopy_indexes_to_remove.size(); i++)
        canopies_to_filter.erase(canopies_to_filter.begin() + canopy_indexes_to_remove[i]);

}

std::vector<Canopy*> CanopyClusteringAlg::multi_core_run_clustering_on(vector<Point*>& points, int num_threads, double max_canopy_dist, double max_close_dist, double max_merge_dist, double min_step_dist, double stop_proportion_of_points, int stop_num_single_point_clusters, string canopy_size_stats_fp, bool show_progress_bar, TimeProfile& time_profile){

    _log(logINFO) << "General:";
    _log(logINFO) << "num_threads:\t " << num_threads;
    _log(logINFO) << "";
    _log(logINFO) << "Algorithm Parameters:";
    _log(logINFO) << "max_canopy_dist:\t " << max_canopy_dist;
    _log(logINFO) << "max_close_dist:\t " << max_close_dist;
    _log(logINFO) << "max_merge_dist:\t " << max_merge_dist;
    _log(logINFO) << "min_step_dist:\t " << min_step_dist;
    _log(logINFO) << "";
    _log(logINFO) << "Early stopping:";
    _log(logINFO) << "stop_proportion_of_points:\t " << stop_proportion_of_points;
    _log(logINFO) << "stop_num_single_point_clusters:\t " << stop_num_single_point_clusters;

    _log(logPROGRESS) << "############ Shuffling ############";
    time_profile.start_timer("Shuffling");
    std::srand ( unsigned ( std::time(NULL) ) );
    std::random_shuffle(points.begin(), points.end());
    time_profile.stop_timer("Shuffling");

    _log(logPROGRESS) << "############ Creating Canopies ############";
    boost::unordered_set<Point*> marked_points;//Points that should not be investigated as origins
    vector<unsigned int> canopy_size_per_origin_num;//Contains size of the canopy created from origin by it's number, so first origin gave canopy of size 5, second origin gave canopy of size 8 and so on
    int num_of_consecutive_canopies_of_size_1 = 0;
    int last_progress_displayed_at_num_points = 0;

    std::vector<Canopy*> canopy_vector;

    vector<Point*> close_points;
    close_points.reserve(points.size());

    



    //
    //Create canopies
    //
    time_profile.start_timer("Clustering");
        
    int num_canopy_jumps = 0;


#pragma omp parallel for shared(marked_points, canopy_vector, num_canopy_jumps, canopy_size_per_origin_num, num_of_consecutive_canopies_of_size_1) firstprivate(close_points, max_canopy_dist, max_close_dist, max_merge_dist, min_step_dist, last_progress_displayed_at_num_points) schedule(dynamic)
    for(int origin_i = 0; origin_i < points.size(); origin_i++){

        //Early stopping proportion of points
        if(marked_points.size() > stop_proportion_of_points * points.size()){
            continue;
        }

        if(num_of_consecutive_canopies_of_size_1 == stop_num_single_point_clusters){
            continue;
        }

        //Show progress bar
        {
            //Only master thread executes this
            if(omp_get_thread_num() == 0){
                if(log_level >= logPROGRESS && show_progress_bar){
                    if(marked_points.size() > last_progress_displayed_at_num_points + stop_proportion_of_points * points.size()/100){
                        printProgBar(marked_points.size(),stop_proportion_of_points * points.size());
                        last_progress_displayed_at_num_points = marked_points.size();
                    }
                }
            }
        }



        Point* origin = points[origin_i]; 

        if(marked_points.find(origin) != marked_points.end())
            continue;

        {
            _log(logDEBUG) << "Unmarked points count: " << points.size() - marked_points.size() << " Marked points count: " << marked_points.size();
            _log(logDEBUG) << "points.size: " << points.size() << " origin_i: " << origin_i << " origin->id: " << origin->id ;

            _log(logDEBUG1) << "Current canopy origin: " << origin->id;
        }

        Canopy *c1;
        Canopy *c2;

        c1 = create_canopy(origin, points, close_points, max_canopy_dist, max_close_dist, true);

        c2 = create_canopy(c1->center, points, close_points, max_canopy_dist, max_close_dist, false);

        double dist = get_distance_between_points(c1->center, c2->center);

        {
            _log(logDEBUG2) << *c1;
            _log(logDEBUG2) << *c2;
            _log(logDEBUG2) << "dist: " << dist;
        }

        {
            _log(logDEBUG3) << "Point1:" ;
            _log(logDEBUG3) << *c1 ;
            _log(logDEBUG3) << "Point2:" ;
            _log(logDEBUG3) << *c2 ;
            _log(logDEBUG3) << "First potential jump correlation: " << dist;
        }
        

        while(dist > min_step_dist){
            delete c1;
            c1=c2;

#pragma omp atomic
            num_canopy_jumps++;

            c2=create_canopy(c1->center, points, close_points, max_canopy_dist, max_close_dist, false);
            dist = get_distance_between_points(c1->center, c2->center); 
            _log(logDEBUG2) << *c1;
            _log(logDEBUG2) << *c2;
            _log(logDEBUG2) << "distance: " << dist;
        }

        //Now we know that c1 and c2 are close enough and we should choose the one that has more neighbours

        Canopy* final_canopy = c1->neighbours.size() > c2->neighbours.size() ? c1 : c2;

#pragma omp critical
        {
            //Do not commit anything if by chance another thread marked the current origin
            if(marked_points.find(origin) == marked_points.end()){

                //Add canopy
                marked_points.insert(origin);

                canopy_vector.push_back(final_canopy);

                BOOST_FOREACH(Point* n, c1->neighbours){
                    marked_points.insert(n);
                }

                //TODO: Could be done better
                if(!(final_canopy->origin->id.compare("!GENERATED!"))){
                    marked_points.insert(c1->origin);
                }

                //Early stopping by number of number of consecutive canopies of size 1
                if(final_canopy->neighbours.size() == 1)
                    num_of_consecutive_canopies_of_size_1++;
                else
                    num_of_consecutive_canopies_of_size_1 = 0;

                if(num_of_consecutive_canopies_of_size_1 >= stop_num_single_point_clusters){
                    _log(logINFO) << "Reached " << num_of_consecutive_canopies_of_size_1  << " of consecutive canopies of size 1. Stopping.";
                }

                //Statistics showing size of canopies per analyzed origin
                canopy_size_per_origin_num.push_back(final_canopy->neighbours.size());

            }

        }

    }
    time_profile.stop_timer("Clustering");

    //Save canopy size statistics if requested
    if(canopy_size_stats_fp != ""){
        try{
            ofstream canopy_size_stats_file;
            canopy_size_stats_file.open(canopy_size_stats_fp.c_str(), ios::out | ios::trunc);
            BOOST_FOREACH(int canopy_size, canopy_size_per_origin_num){
                canopy_size_stats_file << canopy_size << endl;
            }
            canopy_size_stats_file.close();
        } catch (ios_base::failure){
            _log(logWARN) << "Error occured when trying to save canopy statistics to: " << canopy_size_stats_fp << ". I will skip this task.";
        }
    }

    _log(logINFO) << "";
    _log(logINFO) << "Avg. number of canopy jumps: " << num_canopy_jumps/(double)canopy_vector.size();
    _log(logINFO) << "Number of canopies before merging: " << canopy_vector.size();

    int original_number_of_canopies = canopy_vector.size();

    //
    // Merge Canopies
    //
    std::vector<Canopy*> merged_canopy_vector;

    time_profile.start_timer("Merging");
    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "############Merging canopies#############";
    while(canopy_vector.size()){
        bool any_canopy_merged = false;

        std::vector<int> canopies_to_be_merged_index_vector;
        std::vector<Canopy*> canopies_to_merge;

        //This is the canopy we will look for partners for
        Canopy* c = canopy_vector[0];

        //Get indexes of those canopies that are nearby
#pragma omp parallel for shared(canopies_to_be_merged_index_vector) 
        for(int i = 1; i < canopy_vector.size(); i++){


            Canopy* c2 = canopy_vector[i]; 

            _log(logDEBUG2) << "Calculating distances";
            _log(logDEBUG2) << *c->center;
            _log(logDEBUG2) << *c2->center;

            double dist = get_distance_between_points(c->center, c2->center);

            _log(logDEBUG2) << "Distance: " << dist;

            if(dist < max_merge_dist){
#pragma omp critical
                {
                    canopies_to_be_merged_index_vector.push_back(i);
                }
            }

        }

        //Assemble all canopies close enough in one bag 
        if( canopies_to_be_merged_index_vector.size() ){

            //Put all canopies to be merged in one bag
            canopies_to_merge.push_back(canopy_vector[0]);

            BOOST_FOREACH(int i, canopies_to_be_merged_index_vector)
                canopies_to_merge.push_back(canopy_vector[i]);

        }
            

        //Remove canopies that were close enough from the canopy vector
        if( canopies_to_be_merged_index_vector.size() ){
            //Now we delete those canopies that we merged

            //First - sort so that the indexes won't change during our merger
            std::sort(canopies_to_be_merged_index_vector.begin(), canopies_to_be_merged_index_vector.end());

            while(canopies_to_be_merged_index_vector.size()){
                canopy_vector.erase(canopy_vector.begin() + canopies_to_be_merged_index_vector.back());
                canopies_to_be_merged_index_vector.pop_back();
            }
            canopy_vector.erase(canopy_vector.begin());

            
        }

        //Actually merge the canopies and add the new one to the beginning of the canopy vector
        if( canopies_to_merge.size() ){
            any_canopy_merged = true;

            Canopy* merged_canopy= new Canopy();

            BOOST_FOREACH(Canopy* canopy, canopies_to_merge){
               merged_canopy->neighbours.insert(merged_canopy->neighbours.begin(), canopy->neighbours.begin(), canopy->neighbours.end()); 
                //delete those canopies which are merged into one
               delete canopy;
            }

            merged_canopy->center = get_centroid_of_points( merged_canopy->neighbours );

            canopy_vector.insert(canopy_vector.begin(), merged_canopy);

        }

        //If no canopies were merged remove the canopy we compared against the others
        if( !any_canopy_merged ){
            merged_canopy_vector.push_back(canopy_vector[0]);
            canopy_vector.erase(canopy_vector.begin());
        }

        //Show progress bar
        {
            if(log_level >= logPROGRESS){
                printProgBar(original_number_of_canopies - canopy_vector.size(), original_number_of_canopies );
            }
        }
    }
    time_profile.stop_timer("Merging");

    _log(logINFO) << "Number of canopies after merging: " << merged_canopy_vector.size();

    return merged_canopy_vector;

}

