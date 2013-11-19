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
Last edited: 11/18/13
"""""""""""""""""""""""""""""""""""

def main():
    
    #Create new graph, fill it with some values

    nodeTuples = [(0,0), (1, 0), (2,1), (3,1), (4,0)]
    edgeTuples = [(1,4,1), (0,2,-1), (3,1,-1), (2,4,1)]

    network = gn.geneNetwork(nodeTuples, edgeTuples)
    network.printNetwork('graphBefore.jpg')

    network.update(1.0)

    network.printNetwork('graphAfter.jpg')


if __name__ == "__main__":
    main()
