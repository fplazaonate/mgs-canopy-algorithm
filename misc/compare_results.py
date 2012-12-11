#!/opt/local/bin/python2.7

import sys

reference_results = {}

for line in open(sys.argv[1],"r").readlines():
    words = line.strip().split()
    canopy_num = words[0].split(":")[1]
    gene_id = words[1]
    reference_results.setdefault(canopy_num, set()).add(gene_id)

reference_results = list(reference_results.itervalues())

reference_results.sort(key=lambda x: len(x))


my_results = []

for line in open(sys.argv[2],"r").readlines():
    words = line.strip().split(" \t,")
    gene_ids = set(words[1:]) 
    my_results.append(gene_ids)

my_results.sort(key=lambda x: len(x))

ref_canopy_stats = []
for ref_canopy in reference_results:
    ref_canopy_stats.append(len(ref_canopy))

#compare 
my_canopy_stats = []

for canopy in filter(lambda x: len(x), my_results):
    best_match = 0
    for ref_canopy in reference_results:
        match = len(ref_canopy & canopy)/len(canopy)
        if match > best_match:
            best_match = match

    my_canopy_stats.append((len(canopy), best_match))







