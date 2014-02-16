

""""""""""""""""""""""""""""""""""""""""""""""""""""""""
GIRVAN LAB ROTATION PROJECT: DIVERSITY OF GENE NETWORKS

James Komianos, Last edited: 2/12/2014
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SUMMARY OF FILES:

genenetwork.py: This contains a class to create a general gene network.
Included are constructors, a reset function, and functions to calculate
expression diversities, hamming distances, and attractor counts. 

NOTE: There is currently a slight bug in hamming distance code, something
got changed around and I haven't figured it out yet.

genetest.py: This is a main testing module that contains examples on how
to run the code and generate plots for expression diversity, hamming distance, and attractor count. The three tests given are good examples of 
how to run those three functions. Also included was a timing test for benchmarking the gene expression diversity code.

25NodeTest: This contains plots of diversity parameter vs threshold for a 25 node network, with varying numbers of initial conditions per threshold.

