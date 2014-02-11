#!/usr/bin/env python

import inspect
import os
import sys

def load_from_file(filename):
    graph = {'g': {}, 'n': 0, 'm': 0}
    sparsify_threshold = 0.0
    for line in open(filename):
        id1, id2, w = line.split(' ')

        w = float(w)

        if w > sparsify_threshold:

            if not id1 in graph['g']:
                graph['g'][id1] = {}
                graph['n']+=1
            if not id2 in graph['g']:
                graph['g'][id2] = {}
                graph['n']+=1

            graph['g'][id1][id2] = w
            graph['g'][id2][id1] = w

            graph['m']+=1

    return graph

def convert_to_pajek(filename):

    graph = load_from_file(filename)
    output_filename = "/mnt/{}.pajek.net".format(filename)
    out = open(output_filename, "w")
    out.write("*Vertices " + str(graph['n']) + "\n")

    for i in range(graph['n']):
        out.write(str(i+1) + " \"" + str(i) + "\" 0.0 0.0 0.0\n")

    out.write("*Edges\n")

    for i in graph['g']:
        si = int(i)
        for n in graph['g'][i]:
            sj = int(n)
            if si < sj:
                # print only one edge, assumes undirected networks
                out.write(str(si+1) + " " + str(sj+1) + " " + str(graph['g'][str(si)][str(sj)]) + "\n")
    out.close()

    return output_filename

def run_clustering(filename, num_cores = 32):
    pwd = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    ## compile
    os.system("gcc -O2 -Wall -o {}/../clusterPBD {}/../clusterPBD.c -lm -lpthread".format(pwd, pwd))
    ## convert the distance matrix (Wreduced.txt) to pajek format (expected by clusterPBD)
    output_filename = convert_to_pajek(filename)
    ## run the algorithm
    os.system("{}/../clusterPBD {} 0 2 0 {}".format(pwd, output_filename, num_cores))

    ## the results are partition_N_clusters.class for N in [9, 50, 100, 200, 500, 1000, 5000, 10000]

    print "Done."
    

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print ""
        print "Usage: "
        print "         run_clustering name_of_the_distance_matrix [number_of_cores]"
        print ""
    else:
        run_clustering(sys.argv[0])
