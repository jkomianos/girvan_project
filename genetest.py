import math
import time
import sys
import scipy as sp
import numpy as np
sys.path.append("~/work/girvan_project/")
import genenetwork as gn
import matplotlib.pyplot as plt

"""""""""""""""""""""""""""""""""""
Code to test the geneNetwork class.

James Komianos
Last edited: 1/29/14
"""""""""""""""""""""""""""""""""""

def timingGEDTest():

    #Testing expression diversity code
    nStart = 200
    nEnd = 1600

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
    plt.plot(xrange(nStart, nEnd, 100), times)
    plt.ylabel('Time to compute (s)')
    plt.xlabel('Number of nodes') 
    plt.grid(True)
    plt.savefig("GEDBenchmark4.jpg")

def periodicOrbitTest():

    numGenes = 25
    network = gn.geneNetwork(preset="RANDOM", numNodes=numGenes, connectionsPerNode=int(0.2*numGenes))

    #print network.findPeriodicOrbit(threshold = 0)
    #print  network.countAttractors(numICs=100, threshold=0)

    network.generateAttractorsVsThreshold(thresholdMin=-10, thresholdMax=15,
                                        thresholdStep=0.25, numICsPerThreshold=400)


def geneExpressionDiversityTest():

    numGenes = 25
    network = gn.geneNetwork(preset="RANDOM", numNodes=numGenes, connectionsPerNode=int(0.2*numGenes))
    network.generateExpDivVsThreshold(thresholdMin=-10, thresholdMax=15,
                                        thresholdStep=0.25, numICsPerThreshold = 400)

def hammingDistanceTest():
    numGenes = 1000
    network = gn.geneNetwork(preset="RANDOM", numNodes=numGenes, connectionsPerNode=int(0.2*numGenes))
    network.generateHammingDistVsThreshold(thresholdMin=-4, thresholdMax=5,
                                        thresholdStep=0.2, numICsPerThreshold = 10, 
                                        numTimeStepsToSaturate=10, diffPermutation=0.1)

def diversityTest():
    numGenes = 25
    network = gn.geneNetwork(preset="RANDOM", numNodes=numGenes, connectionsPerNode=int(0.2*numGenes))
    network.generateDiversityVsThreshold(thresholdMin=-3, thresholdMax=4,
                                        thresholdStep=0.1, numICsPerThreshold = 200)


def main():
    
    #timingGEDTest()
    #periodicOrbitTest()
    #geneExpressionDiversityTest()
    #hammingDistanceTest()
    diversityTest()


if __name__ == "__main__":
    main()
