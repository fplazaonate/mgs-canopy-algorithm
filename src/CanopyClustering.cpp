#include <algorithm>


#include <boost/foreach.hpp>

#include <CanopyClustering.hpp>
#include <Log.hpp>

Canopy* CanopyClusteringAlg::create_canopy(Point* origin, boost::unordered_set<Point*>& marked_points, std::vector<Point*>& points, double min_correlation){

    std::vector<Point*> neighbours;

    //TODO: dangerous?
    BOOST_FOREACH(Point* potential_neighbour, points){
        if(Point::get_distance_between_points(origin, potential_neighbour) > min_correlation){
            neighbours.push_back(potential_neighbour);
        }
    }

    Point* center;

    if(!neighbours.size())
        center = origin;
    else
        center = Point::get_centroid_of_points(neighbours);

    return new Canopy(origin, center, neighbours);

}

std::vector<Canopy*> CanopyClusteringAlg::single_core_run_clustering_on(std::vector<Point*>& points){

    _log(logDEBUG1) << "############Creating Canopies############";

    double min_canopy_correlation = 0.9;
    double canopy_merge_distance_threshold = 0.97;
    double canopy_iteration_min_correlation = 0.1;

    //TODO: ensure it is actually random
    std::random_shuffle(points.begin(), points.end());

    //TODO: this will be super slow!
    boost::unordered_set<Point*> marked_points;

    std::vector<Canopy*> canopy_vector;

    //
    //Create canopies
    //
    BOOST_FOREACH(Point* origin, points){
        
        _log(logINFO) << "Unmarked points count: " << points.size() - marked_points.size();

        if(marked_points.find(origin) != marked_points.end())
            continue;

        _log(logDEBUG2) << "Current canopy origin: " << origin->id;

        Canopy *c1;
        Canopy *c2;

        c1 = create_canopy(origin, marked_points, points, min_canopy_correlation );
        if(c1->neighbours.size())
            c2->origin = Point::get_centroid_of_points(c1->neighbours);
        else
            c2->origin = c1->origin;

        c2 = create_canopy(c2->origin, marked_points, points, min_canopy_correlation);

        double correlation = Point::get_distance_between_points(c1->center, c2->center);

        _log(logDEBUG3) << *c1;
        _log(logDEBUG3) << *c2;
        _log(logDEBUG3) << "correlation: " << correlation;


        while(correlation < canopy_iteration_min_correlation){
            c1=c2;

            c2=create_canopy(c1->center, marked_points, points, min_canopy_correlation);
            correlation = Point::get_distance_between_points(c1->center, c2->center); 
            _log(logDEBUG3) << *c1;
            _log(logDEBUG3) << *c2;
            _log(logDEBUG3) << "correlation: " << correlation;
        }

        canopy_vector.push_back(c2);

        BOOST_FOREACH(Point* n, c2->neighbours){
            marked_points.insert(n);
        }

    }

    std::vector<Canopy*> merged_canopy_vector;

    //Merge canopies
    //
    //The basic idea:
    //1. We pick one canopy c from the canopy vector (canopy with index 0)
    //2. We go through remaining canopies in that vector
    //3. If any of the remaining canopies is closer to c than threshold
    //4. We 
    
    _log(logDEBUG1) << "############Merging canopies############";
    while(canopy_vector.size()){
        bool any_canopy_merged = false;

        std::vector<int> canopies_to_be_merged_index_vector;
        std::vector<Canopy*> canopies_to_merge;

        //This is the canopy we will look for partners for
        Canopy* c = canopy_vector[0];

        //Get indexes of those canopies that are nearby
        for(int i = 1; i < canopy_vector.size(); i++){

            Canopy* c2 = canopy_vector[i]; 

            _log(logDEBUG3) << "Calculating distances";
            _log(logDEBUG3) << *c->center;
            _log(logDEBUG3) << *c2->center;

            double correlation = Point::get_distance_between_points(c->center, c2->center);

            _log(logDEBUG3) << "Correlation: " << correlation;

            if(correlation > canopy_merge_distance_threshold)
                canopies_to_be_merged_index_vector.push_back(i);

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

            Canopy* merged_canopy;

            BOOST_FOREACH(Canopy* canopy, canopies_to_merge)
               merged_canopy->neighbours.insert(merged_canopy->neighbours.begin(), canopy->neighbours.begin(), canopy->neighbours.end()); 

            merged_canopy->center = Point::get_centroid_of_points( merged_canopy->neighbours );

            canopy_vector.insert(canopy_vector.begin(), merged_canopy);
        }

        //If no canopies were merged remove the canopy we compared against the others
        if( !any_canopy_merged ){
            merged_canopy_vector.push_back(canopy_vector[0]);
            canopy_vector.erase(canopy_vector.begin());
        }
    }
    return merged_canopy_vector;

    //return canopy_vector;

}
