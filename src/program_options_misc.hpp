#include <iostream>
#include <fstream>
#include <limits>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/join.hpp>

using namespace std;

bool check_if_within_bounds(string option_name, double value, double lower, double higher){
    if( value >= lower - numeric_limits<double>::epsilon() && value <= higher + numeric_limits<double>::epsilon()){
        return true;
    }else{ 
        _log(logERR) << "Option: \"" << option_name << "\" must be a value within range: <" << lower << ";" << higher << ">";
        exit(1);
    }
}

bool check_if_within_bounds(string option_name, int value, int lower, int higher){
    if( value >= lower && value <= higher ){
        return true;
    }else{ 
        _log(logERR) << "Option: \"" << option_name << "\" must be a value within range: <" << lower << ";" << higher << ">";
        exit(1);
    }
}

bool check_if_file_is_readable(string option_name, string path){
    ofstream file;
    try{
        file.open(path.c_str(), ios::in );
        file.close();
        return true;
    } catch (ios_base::failure){
        _log(logERR) << "Option: \"" << option_name << "\" must be accessible and readable.";
        exit(1);
    }
}

bool check_if_file_is_writable(string option_name, string path){
    ofstream file;
    try{
        file.open(path.c_str(), ios::out | ios::trunc);
        file.close();
        return true;
    } catch (ios_base::failure){
        _log(logERR) << "Option: \"" << option_name << "\" must be accessible and writable.";
        exit(1);
    }
}

bool check_if_one_of(string option_name, string value, vector<string> valid_options){
    BOOST_FOREACH(string valid_opt, valid_options){
        if( value == valid_opt){
            return true;
        }
    }
    _log(logERR) << "Option: \"" << option_name << "\" must be one of:" << boost::algorithm::join(valid_options, ", ");
    exit(1);
}
    

