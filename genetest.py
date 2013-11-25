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

    network = gn.geneNetwork(preset="RANDOM")

    network.printNetwork('graph0.jpg')
    network.testUpdate()
    network.printNetwork('graph1.jpg')
    network.testUpdate()
    network.printNetwork('graph2.jpg')
    network.testUpdate()
    network.printNetwork('graph3.jpg')
    network.testUpdate()
    network.printNetwork('graph4.jpg')



if __name__ == "__main__":
    main()
