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
Last edited: 1/13/14
"""""""""""""""""""""""""""""""""""

def main():
    
    #Testing expression diversity code
    nStart = 200
    nEnd = 1000

    times = []

    for numGenes in xrange(nStart, nEnd, 100):

        network = gn.geneNetwork(preset="RANDOM", numNodes=numGenes, connectionsPerNode=int(0.01*numGenes))

        t1 = time.time()
        network.generateExpDivVsThreshold(thresholdMin=0, thresholdMax=1,
                                            thresholdStep=1, numICsPerThreshold=10)
        t2 = time.time()

        times.append(t2 - t1)

    #Plot the timings
    plt.clf()
    plt.plot(times, xrange(nStart, nEnd, 100))
    plt.xlabel('Time to compute')
    plt.ylabel('Number of nodes') 
    plt.grid(True)
    plt.savefig("GEDBenchmark.jpg")


if __name__ == "__main__":
    main()
