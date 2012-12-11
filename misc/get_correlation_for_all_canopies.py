import sys
import pickle
from scipy import array
from scipy.stats.stats import pearsonr
from random import choice

reference_results_file = "/Users/piotr/Documents/projects/canopy_clustering/data/C2GU_sets.txt"
points_file = "/Users/piotr/Documents/projects/canopy_clustering/data/MetaHIT663.RefGeneCat396.basecov_v2"

#Import reference results
reference_results = {}

for line in open(reference_results_file,"r").readlines():
    words = line.strip().split()
    canopy_num = words[0].split(":")[1]
    gene_id = words[1]
    reference_results.setdefault(canopy_num, set()).add(gene_id)
    
reference_results_arr = []
for genes_set in reference_results.itervalues():
    reference_results_arr.append(frozenset(genes_set))

genes_from_some_random_bjorn_canopies = set()


#Read point sample data
point_sample_data = {}

for line in open(points_file,"r"):
    words = line.strip().split()
    if words[0] in genes_from_some_random_bjorn_canopies:
        point_sample_data[words[0]] = array(map(lambda x: float(x), words[1:]))

correlations = []
for canopy_num, canopy in enumerate(reference_results_arr):
    
    print str(canopy_num) + ": current canopy length: " + str(len(canopy))

    #Calculate Centroid
    centroid_list = []
    genes = list(canopy)
    num_data_samples = len(point_sample_data[genes[0]])
    
    for i in range(num_data_samples):
        column = []
        for gene in canopy:
            column.append(point_sample_data[gene][i])
        column.sort()
        mean = column[int(floor(len(canopy)/2))]
        centroid_list.append(mean)
        
    centroid_vec = array(centroid_list)
    
    #Calculate correlations between points and centroid
    for gene in canopy:
        gene_vec = point_sample_data[gene]
        correlation = pearsonr(gene_vec,centroid_vec)[0]
        
        correlations.append((gene, correlation, centroid_vec, gene_vec))

#Show histogram of correlations
correlation_values = map(lambda x: x[1], correlations)

pickle.dump(correlation_values, "get_correlation_for_all_canopies.out.pickle")


