import math
import time
import sys
import scipy as sp
import numpy as np
sys.path.append("~/work/girvan_project/")
import genenetwork as gn

"""""""""""""""""""""""""""""""""""
Code to test the geneNetwork class.

James Komianos
Last edited: 12/12/13
"""""""""""""""""""""""""""""""""""

def main():
    
    #Create new graph, fill it with some values

    network = gn.geneNetwork(preset="RANDOM", numNodes=1000, connectionsPerNode=10)

    network.generateHammingDistVsThreshold(thresholdMin=-2, thresholdMax=6,
                                        thresholdStep=0.25, numICsPerThreshold=100,
                                        numTimeStepsToSaturate=10, diffPermutation = 0.01)



if __name__ == "__main__":
    main()
