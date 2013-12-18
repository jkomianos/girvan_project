import scipy as sp
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import math
import time
import sys
import copy

"""""""""""""""""""""""""""""""""""""""""""""""
A geneNetwork class for gene network analysis.

James Komianos
Last Edited: 12/18/13
"""""""""""""""""""""""""""""""""""""""""""""""

class geneNetwork:

	G = nx.Graph()
	GPerm = None

	#To keep track of node values
	GNodes = {}
	GPermNodes = {}

	numGenes = 0
	connectionsPerGene = 0

	"""
	Constructor
	nodes= a list tuples containing node name and value
	edgeTuples= a list of tuples containing
		the edges to add, with weights
	"""
	def __init__(self, preset="RANDOM", numNodes=1000, 
				connectionsPerNode=10, nodeTuples=None, edgeTuples=None):

		if preset == "RANDOM":
			#Add nodes
			for i in xrange(0, numNodes):

				state = 1 if random.random() >= 0.5 else 0
				self.G.add_node(i, value=state)
				self.GNodes[i] = state

			#Now add edges
			for i in xrange(0, numNodes):

				#Check how many connections this node already has
				edgesToAdd = connectionsPerNode - len(self.G.neighbors(i))

				for j in xrange(0, edgesToAdd):
					randomNode = random.randint(0,numNodes-1)
					randomWeight = 1 if random.random() >= 0.33 else -1
					toConnect = randomNode if randomNode != i else random.randint(0,numNodes-1)
					self.G.add_edge(i, randomNode, weight=randomWeight)

			self.numGenes = numNodes
			self.connectionsPerGene = connectionsPerNode

		if preset == "MANUAL":
			for nodeTuple in nodeTuples:
				self.G.add_node(nodeTuple[0], value=nodeTuple[1])
				self.GNodes[nodeTuple[0]] = nodeTuple[1]

			self.G.add_weighted_edges_from(edgeTuples)
			self.numGenes = len(nodeTuples)


	"""
	Resets the current graph	
	"""
	def reset(self, resetEdges=True):

		if resetEdges:

			self.G = nx.Graph()
			self.GPerm = None

			#Add nodes
			for i in xrange(0, self.numGenes):
				state = 1 if random.random() >= 0.5 else 0
				self.G.add_node(i, value=state)
				self.GNodes[i] = state

			#Now add edges
			for i in xrange(0, self.numGenes):
				for j in xrange(0, self.connectionsPerGene):
					randomNode = random.randint(0,self.numGenes-1)
					randomWeight = 1 if random.random() >= 0.33 else -1
					toConnect = randomNode if randomNode != i else random.randint(0,self.numGenes-1)
					self.G.add_edge(i, randomNode, weight=randomWeight)

		#Just reset node values			
		else:
			for n in self.G:
				state = 1 if random.random() >= 0.5 else 0
				self.G.node[n]['value']	= self.GNodes[n] = state

	"""
	Update graphs using thresholding
	"""	
	def update(self, threshold):

		GNodesChanged = []
		GPermNodesChanged = []

		#Iterate over all nodes
		for n in self.G:
			sumN = 0
			#Find all adjacents
			for neighbor in self.G.neighbors_iter(n):
				#Get edge weight
				value = self.GNodes[neighbor]
				weight = self.G[n][neighbor]['weight']
				sumN += value * weight	

			#Update node based on sign of sum
			if (sumN > threshold): 
				self.G.node[n]['value'] = 1
				if self.GNodes[n] != 1: GNodesChanged.append(n)
			else: 
				self.G.node[n]['value'] = 0
				if self.GNodes[n] != 0: GNodesChanged.append(n)

		#Update the permutation if it exists
		if self.GPerm != None:	
			for n in self.GPerm:
				sumN = 0
				#Find all adjacents
				for neighbor in self.GPerm.neighbors_iter(n):
					#Get edge weight 
					value = self.GPermNodes[neighbor]
					weight = self.GPerm[n][neighbor]['weight']
					sumN += value * weight

				#Update node based on sign of sum	
				if (sumN > threshold): 
					self.GPerm.node[n]['value'] = 1 
					if self.GPermNodes[n] != 1: GNodesChanged.append(n)
				else: 
					self.GPerm.node[n]['value'] = 0
					if self.GPermNodes[n] != 0: GPermNodesChanged.append(n)

		#Update the changed nodes in the dictionary
		for n in GNodesChanged:
			currentVal = self.GNodes[n]
			self.GNodes[n] = 1 if currentVal == 0 else 0	

		for n in GPermNodesChanged:
			currentVal = self.GPermNodes[n]
			self.GPermNodes[n] = 1 if currentVal == 0 else 0
					
	"""
	Create permutation of current graph
	"""
	def createPermutation(self, diffThreshold=0.01):
		self.GPerm = self.G.copy()
		#Change percentage of nodes (keep edges the same)
		numToChange = int(diffThreshold * self.numGenes)
		self.GPermNodes = copy.deepcopy(self.GNodes)

		for i in xrange(0, numToChange):
			nodeToChange = random.randint(0,self.numGenes - 1)
			currentVal = self.GPerm.node[nodeToChange]['value']
			state = 1 if currentVal == 0 else 0
			self.GPerm.node[nodeToChange]['value'] = self.GPermNodes[nodeToChange] = state


	"""
	Compute the hamming distance between two graphs (G and GPerm)
	"""		
	def calculateHammingDist(self, other=None):

		hammingDist = 0

		if other == None:
			for i in xrange(0, self.numGenes):
				if self.G.node[i]['value'] != self.GPerm.node[i]['value']:
					hammingDist += 1
		else:
			for i in xrange(0, self.numGenes):
				if self.G.node[i]['value'] != other.node[i]['value']:
					hammingDist += 1
		
		return hammingDist

	"""
	Compute the saturated hamming distance between G and GPerm
	"""
	def calculateSatHammingDist(self, numTimeSteps=10, threshold=0):	

		for t in xrange(0, numTimeSteps):
			#Update, calculate hamming distance
			self.update(threshold)

		return float(self.calculateHammingDist()) / self.numGenes

	"""
	Compute a set of data representing avg saturated hamming 
	dist vs thresholds, plot it
	"""
	def generateHammingDistVsThreshold(self, thresholdMin=0, thresholdMax=10,
										thresholdStep=1, numICsPerThreshold=10,
										numTimeStepsToSaturate=10, diffPermutation = 0.01):
	
		self.createPermutation(diffThreshold=diffPermutation)

		thetas = []
		avgHammingDists = []
		for t in np.arange(thresholdMin, thresholdMax, thresholdStep):

			hammingDistsForT = []
			for i in xrange(0, numTimeStepsToSaturate):
				#Create permutation, calculate hammingDist after 
				hammingDistsForT.append(self.calculateSatHammingDist(numTimeSteps=numTimeStepsToSaturate, threshold=t))	
	    		self.reset(resetEdges=True)	
	    		self.createPermutation(diffThreshold=diffPermutation)

			avgHammingDists.append(sum(hammingDistsForT) / float(len(hammingDistsForT)))
			thetas.append(t)

	    #Now, plot this using matplotlib
		plt.clf()
		plt.plot(thetas, avgHammingDists)
		plt.xlabel('Threshold')
		plt.ylabel('Saturated hamming distance (normalized)')	
		plt.grid(True)
		plt.savefig("satHammingDist.jpg")

	"""
	Create gene expression factor for a given network
	"""	
	def createGeneExpFactor(self, threshold):

		#Starting and stopping steps for construction
		start = int(0.5 * self.numGenes)
		stop = int(2 * self.numGenes)

		X = [0] * self.numGenes #gene expression factor

		for t in xrange(0, stop):
			#Do not add to expression factor
			if t < start:
				self.update(threshold=threshold)
			#Add to expression factor	
			else:
				self.update(threshold=threshold)
				for n, value in self.GNodes.iteritems():
					X[n] += value

		#Now, normalize
		return [float(value) / (stop - start) for value in X]					
 
	"""
	Calculate expression diversity for this network	
	"""
	def expressionDiversity(self, threshold=0, numICs=100):
	
		allX = []
		expDiv = 0
		#Create expression factors
		for i in xrange(0, numICs):
			allX.append(self.createGeneExpFactor(threshold=threshold))
			self.reset(resetEdges=False)
			print allX[i]

		#Now, sum over all expression factors
		for m in xrange(0, numICs):
			#Loop over all others
			for n in xrange(0, numICs):
				if m != n:
					for j in xrange(0, self.numGenes): 
						expDiv = allX[m][j] - allX[n][j]

		#normalize and return
		return expDiv / (self.numGenes * numICs)


	"""
	For debugging purposes
	"""

	def testUpdate(self, threshold=0):
		gOriginal = self.G.copy()

		self.update(threshold)
		print self.calculateHammingDist(other=gOriginal)

	def printNetwork(self, outputFile):

		plt.clf()
		pos=nx.spring_layout(self.G)

		val_map = {1: 1.0, 0: 0.2654}
		values = [val_map.get(self.G.node[n]['value'], 0.1) for n in self.G]

		nx.draw(self.G, pos, cmap=plt.get_cmap('prism'), node_color=values)
		plt.savefig(outputFile)

	def printState(self):
		print "PRINTING CURRENT NODE DATA"
		for i in xrange(0, self.numGenes):
			print "Node=", i, "Value G=", self.G.node[i]['value'], "Value GPerm=", self.GPerm.node[i]['value']




		