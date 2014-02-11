
# ClusterPBD

This is the implementation of the algorithm ClusterPBD published in Physical Review E

*Reference:* Pujol, J.M., Béjar, J. and Delgado, J. "Clustering Algorithm for Determining Community Structure in Large Networks". 
Physical Review E 74 (2007):016107 [Paper](http://cscs.umich.edu/%7Ejmpujol/public/papers/spectralCluster.pdf)

#Build from source

The source code is self-contained, no external dependencies. You can build the binary using gcc:

    gcc -O2 -Wall -o clusterPBD clusterPBD.c -lm

This will create the executable clusterPBD

#Datasets

In data/ you can find several networks

* zachary.net : the typical Karate Club (very small, 34 nodes)
* erdos02b.net : scientist that collaborated with Erdos (small, 6927 nodes)
* condmat.net : social network from physicists (medium, 27519 nodes)

*ClusterPBD* can analyze networks of 500K nodes less than one hour and medium networks like condmat in few seconds.

#Usage

Several examples (fromt the networks in data/)

    ./clusterPBD data/zachary.net 0 2 0
    ./clusterPBD data/erdos02b.net 1 2 0
    ./clusterPBD data/condmat.net 0 2 0

##Parameters

You can get them from the binary if you execute it without parameters.


    ./clusterPBD filename type distance msave

      filename (network in pajek format)
      type {0,1} (pajek format: 0, pajek compact format: 1)
      distance (for the initial seed)
      msave {0,1} (don't save intermediate results: 0, save them: 1)


##Output

The executable will generate 4 output files (in the same directory thant the binary, it's hardcoded :/)

*results.class*

Partition with the modularity value, the last partition is the one you probably one, after the
last "----", each line is the class of the node i-th, finally, the modularity before, modularity
after and finally the number of classes in the partition. 

*results.dendro*

The dendrogram (usable in matlab).

*results.mod*

Modularity by iteration of the PBD algorithm.

*results.join*

The joins of the agglomerative algorithm (after the initial seed)



