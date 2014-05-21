
########################
# Table of contents
########################
1. Introduction
2. Obtaining the code 
3. Obtaining the binary 
4. Run Requirements
5. Input 
6. Output
7. Program parameters and example run command

########################
# 1. Introduction
########################
Our variation of canopy clustering focuses on efficient clustering of points in multi-dimensional pearson correlation space. The basic notion of the heuristic is choosing a point(seed point) at random, and upon establishing it's distance to all the other points, single out those that are within a specified canopy distance. A median profile of those points is then calculated creating a canopy centroid. The canopy creation process is then repeated perpetually(canopy walk) from the previously created centroid until the distance between previous centroid is small enough. All points from the last canopy are then marked, and will not become seed points again. The above will be repeated for all possible seed points.

########################
# 2. Obtaining the code 
########################
It is suggested that the latest copy of the software is obtained from the public code repository at: http://git.dworzynski.eu/mgs-canopy-algorithm/wiki/Home. An example command for doing so is:

$ git clone https://bitbucket.org/HeyHo/mgs-canopy-algorithm.git 

########################
# 3. Obtaining the binary 
########################
Software requirements:

The Canopy Clustering program is implemented entirely in C++. It's only dependencies is BOOST headers and BOOST program options library.
On debian based systems the package names should be libboost-dev and libboost-program-options-dev.

Compiling:
Two Makefiles are provided with the code, for typical x86_64 system please use Makefile-cges2 in the src directory of the downloaded code.

The easiest way to compile the code is to change your current directory to src subdirectory run:

$ make -f Makefile-cges2

If the compilation is successfull a cc.bin executable should appear in the src subdirectory.

########################
# 4. Run Requirements
########################
Our program is suitable for running only on shared-memory machines. By default the dependencies are compiled statically into the binary, but in case your executable was compiled dynamically you must ensure that boost program options library is present on your system and accessible from LD_LIBRARY_PATH or equivalent.

RAM Memory
When run with a matrix of ~7M points with 600 data samples each using 12 threads, the program may take up to 40GB of active memory. The memory requirement scales linearly with both number of points, number of samples and number of threads.

CPU
The program was written with parallelization as it's primary focus and it's running time should scale linearly even with high amounts of cores. There is no real minimum to the number of threads that the program should be run on but we recommend using at least 8 logical-threads(as returned by your /etc/cpuinfo).

Disk
The typical run will generate no more than one-tenth of the input size as it's output.

########################
# 5. Input 
########################

The program expects as an input a file containing points and their profiles.

Input file requirements:
 > points must be specified in a line-by-line fashion
 > the first column must always be point's name
 > second and following columns must be point's dimension-positions (profiles)
 > all columns are whitespace delimited (tabs or spaces)
 > input file cannot have any kind of header all points must have equal amount of data points (same profile lengths)

Example:
Point0 0 0 0 1 3 3
Point1 1 3 0 1 9 6
Point2 0 6 0 1 0 3
Point3 0 3 0 2 2 1
Point4 0 2 0 1 1 7

########################
# 6. Output
########################

The output is composed of two files: the cluster file and the cluster profiles file

Cluster file:
The cluster file contains, line by line, tab separated pairs of values <cluster name> <point name>. Clusters are sorted according to their size from the biggest to the smallest. Points are not sorted in any particular order.

Cluster file example:
MGU00000        Point0
MGU00000        Point3
MGU00000        Point4
MGU00000        Point5
MGU00000        Point7
MGU00001        Point0
MGU00001        Point3
MGU00001        Point8
MGU00001        Point9
MGU00002        Point13
MGU00002        Point12
MGU00002        Point3
MGU00003        Point1
MGU00003        Point5
MGU00004        Point5
MGU00004        Point6
MGU00005        Point8

Cluster profiles file:
Cluster profiles file contains line-by-line the cluster name and it's profile. The first column is always the name and following columns are the profile data-points in the order corresponding to point input file. All columns are tab separated.

Cluster profiles file example:
MGU00000        1 0 0 3 3 3
MGU00001        0 8 0 1 2 1
MGU00002        0 3 0 1 3 3
MGU00003        0 2 1 1 6 3
MGU00004        9 0 0 1 3 0
MGU00005        8 0 0 1 3 0

########################
# 7. Program parameters and example run command
########################
An up-to-date list of program parameters can be obtained by running it with '--help' option which will provide a listing of all available options and their descriptions.

An example of successful run command is presented below. Arguments in "<>" brackets are absolute paths to files.

$ cc.bin -n 16 -i <matrix_input_file> -o <clusters_out> -c <profiles_out> -p <CAG_prefix> --max_canopy_dist 0.1 --max_close_dist 0.4 --max_merge_dist 0.1 --min_step_dist 0.005 --stop_fraction 1 --canopy_size_stats_file <progress_stat_file>

The process can be interrupted using an interrupt system signal. By default the interrupt signal will cause the program to stop the clustering and commence canopy merging - thus the partial results will be provided. This behavior can be changed through "die_on_kill" option.

An example of sending the interrupt signal (SUSE linux):

$ kill -INT [cc.bin PID]






