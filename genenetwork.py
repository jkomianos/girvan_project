import scipy as sp
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import math
import time
import sys
import copy
import itertools
import string


"""""""""""""""""""""""""""""""""""""""""""""""
A geneNetwork class for gene network analysis.

James Komianos
Last Edited: 1/29/14
"""""""""""""""""""""""""""""""""""""""""""""""

class geneNetwork:

	G = None
	GPerm = None

	#To keep track of node values for updating
	#Hold[value, Frozen(T/F)]
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

		self.G = nx.Graph()

		if preset == "RANDOM":
			#Add nodes
			for i in xrange(0, numNodes):

				state = 1 if random.random() >= 0.5 else 0
				self.G.add_node(i, value=state)
				self.GNodes[i] = [state, False]

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
				self.GNodes[nodeTuple[0]] = [nodeTuple[1], False]

			self.G.add_weighted_edges_from(edgeTuples)
			self.numGenes = len(nodeTuples)


	"""
	Finds frozen nodes in graph, removes them from dict
	"""
	def findFrozenNodes(self, threshold):

		for n in self.G:

			weights = []	
			for neighbor in self.G.neighbors_iter(n):
				weights.append(self.G[n][neighbor]['weight'])

			numEdges = len(weights)	
			frozen = True
			truthTable = list(itertools.product([0,1], repeat = numEdges))

			#Determine if node is frozen
			for values in truthTable:
				sumN = 0
				for i in xrange(0, numEdges):
					sumN += values[i] * weights[i]

				if (sumN > threshold):
					frozen = False
					break

			#If in fact frozen, remove from GNodes
			if frozen: self.GNodes[n][1] = True	

		#Remove from permutation if it exists	
		if self.GPerm != None:

			for n in self.GPerm:

				weights = []	
				for neighbor in self.GPerm.neighbors_iter(n):
					weights.append(self.GPerm[n][neighbor]['weight'])

				numEdges = len(weights)	
				frozen = True
				truthTable = list(itertools.product([0,1], repeat = numEdges))

				#Determine if node is frozen
				for values in truthTable:
					sumN = 0
					for i in xrange(0, numEdges):
						sumN += values[i] * weights[i]

					if (sumN > threshold):
						frozen = False
						break

				#If in fact frozen, remove from GNodes
				if frozen: self.GPermNodes[n][1] = True		


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
				self.GNodes[i][0] = state

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
				self.G.node[n]['value']	= self.GNodes[n][0] = state	

	"""
	Update graphs using thresholding
	"""	
	def update(self, threshold):

		#Iterate over all nodes
		for n in self.GNodes.keys():

			if (self.GNodes[n][1] != True):
				sumN = 0
				#Find all adjacents
				for neighbor in self.G.neighbors_iter(n):
					#Get edge weight
					value = self.GNodes[neighbor][0]
					weight = self.G[n][neighbor]['weight']
					sumN += value * weight	

				#Update node based on sign of sum
				valueN = self.GNodes[n][0]
				if (sumN > threshold): 
					self.G.node[n]['value'] = 1
					if valueN != 1: self.GNodes[n][0] = 1
				else: 
					self.G.node[n]['value'] = 0
					if valueN != 0: self.GNodes[n][0] = 0

		#Update the permutation if it exists
		if self.GPerm != None:	
			for n in self.GPermNodes.keys():

				if (self.GPermNodes[n][1] != True):

					sumN = 0
					#Find all adjacents
					for neighbor in self.GPerm.neighbors_iter(n):
						#Get edge weight 
						value = self.GPermNodes[neighbor][0]
						weight = self.GPerm[n][neighbor]['weight']
						sumN += value * weight

					#Update node based on sign of sum	
					valuePermN = self.GPermNodes[n][0]
					if (sumN > threshold): 
						self.GPerm.node[n]['value'] = 1 
						if valuePermN != 1: self.GPermNodes[n][0] = 1
					else: 
						self.GPerm.node[n]['value'] = 0
						if valuePermN != 0: self.GPermNodes[n][0] = 0
					
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
			self.GPerm.node[nodeToChange]['value'] = self.GPermNodes[nodeToChange][0] = state

	"""
	Compute the hamming distance between two graphs (G and GPerm)
	"""		
	def calculateHammingDist(self, other=None):

		hammingDist = 0

		if other == None:
			for i in xrange(0, self.numGenes):
				if self.GNodes[i][0] != self.GPermNodes[i][0]:
					hammingDist += 1
		else:
			for i in xrange(0, self.numGenes):
				if self.GNodes[i][0] != other.node[i]['value']:
					hammingDist += 1
		
		return hammingDist

	"""
	Compute the saturated hamming distance between G and GPerm
	"""
	def calculateSatHammingDist(self, numTimeSteps, threshold):	

		#self.findFrozenNodes(threshold=threshold)

		for t in xrange(0, numTimeSteps):
			#Update, calculate hamming distance
			self.update(threshold)

		return float(self.calculateHammingDist()) / self.numGenes

	"""
	Compute a set of data representing avg saturated hamming 
	dist vs thresholds, plot it
	"""
	def generateHammingDistVsThreshold(self, thresholdMin, thresholdMax,
										thresholdStep, numICsPerThreshold,
										numTimeStepsToSaturate, diffPermutation):
	
		self.createPermutation(diffThreshold=diffPermutation)

		thetas = []
		avgHammingDists = []
		for t in np.arange(thresholdMin, thresholdMax, thresholdStep):

			hammingDistsForT = []
			for i in xrange(0, numTimeStepsToSaturate):
				#Create permutation, calculate hammingDist after 
				hammingDistsForT.append(self.calculateSatHammingDist(numTimeSteps=numTimeStepsToSaturate, threshold=t))	
	    		self.reset(resetEdges=False)	
	    		self.createPermutation(diffThreshold=diffPermutation)

			avgHammingDists.append(sum(hammingDistsForT) / float(len(hammingDistsForT)))
			thetas.append(t)

	    #Now, plot this using matplotlib
		plt.clf()
		plt.plot(thetas, avgHammingDists)
		plt.xlabel('Threshold')
		plt.ylabel('Saturated hamming distance (normalized)')	
		plt.grid(True)
		plt.savefig("hammingDist400ICs.jpg")

	"""
	Create gene expression factor for a given network
	"""	
	def createGeneExpFactor(self, threshold):

		#Starting and stopping steps for construction
		start = int(0.5 * self.numGenes)
		stop = int(2 * self.numGenes)

		X = [0] * self.numGenes #gene expression factor

		#self.findFrozenNodes(threshold=threshold)

		for t in xrange(0, stop):
			#Do not add to expression factor
			if t < start: 
				self.update(threshold=threshold)
				
			#Add to expression factor	
			else:
				self.update(threshold=threshold)
				for n, value in self.GNodes.iteritems():
					X[n] += value[0]

		#Now, normalize
		geneExpFactor = [float(value) / (stop - start) for value in X]
		return geneExpFactor			
 
	"""
	Calculate expression diversity for this network	
	"""
	def expressionDiversity(self, threshold, numICs):
	
		allX = []
		expDiv = 0
		#Create expression factors
		for i in xrange(0, numICs):
			allX.append(self.createGeneExpFactor(threshold=threshold))
			self.reset(resetEdges=False)

		#Now, sum over all expression factors
		for m in xrange(0, numICs):
			#Loop over all others
			for n in xrange(0, numICs):
				if m != n:
					for j in xrange(0, self.numGenes): 
						expDiv += abs(allX[m][j] - allX[n][j])

		#normalize and return
		normExpDiv =  expDiv / (self.numGenes * numICs)

		return normExpDiv

	"""
	Compute a set of data representing gene expression 
	diversity vs thresholds, plot it
	"""
	def generateExpDivVsThreshold(self, thresholdMin, thresholdMax,
										thresholdStep, numICsPerThreshold):
	
		thetas = []
		expressionDiversities = []

		for t in np.arange(thresholdMin, thresholdMax, thresholdStep):

			thetas.append(t)
			expressionDiversities.append(self.expressionDiversity(
								threshold=t, numICs=numICsPerThreshold))	
			self.reset(resetEdges=False)

	    #Now, plot this using matplotlib
		plt.clf()
		plt.plot(thetas, expressionDiversities)
		plt.xlabel('Threshold')
		plt.ylabel('Gene Expression Diversity')	
		plt.grid(True)
		plt.savefig("GED400ICs.jpg")


	"""
	Produces a bitstring of the current state
	"""
	def getStateBitstring(self):	
		return ''.join([str(val[0]) for val in self.GNodes.values()])

	"""
	Given an initial condition, find a periodic orbit and save it
	"""
	def findPeriodicOrbit(self, threshold):

		#Hash table to store the orbit
		orbitTable = {}
		#Find frozen nodes
		#self.findFrozenNodes(threshold = threshold)

		orbitFound = False
		start = True
		currentState = ""
		nextState = ""
		while not orbitFound:

			currentState = self.getStateBitstring() if start else nextState
			self.update(threshold = threshold)
			nextState = self.getStateBitstring()

			#Check if nextState is in orbit
			if nextState in orbitTable:
				orbitFound = True
			else:
				orbitTable[currentState] = nextState
			start = False

		#We have found the orbit, now save it (starting at nextState)
		finalOrbit = []
		saveOrbit = True
		orbitState = firstOrbitState = nextState

		#Save the orbit
		while saveOrbit:
			finalOrbit.append(orbitState)
			try:
				orbitState = orbitTable[orbitState]
				if orbitState == firstOrbitState: saveOrbit = False
			except KeyError:
				saveOrbit = False

		return finalOrbit

	"""
	Counting attractors for a given number of ICs
	"""
	def countAttractors(self, numICs, threshold):

		attractors = []
		for i in xrange(0, numICs):
			orbit = self.findPeriodicOrbit(threshold=threshold)
			if orbit not in attractors:
				attractors.append(orbit)
			self.reset(resetEdges=False)

		return len(attractors)

	"""
	Create attractor count vs threshold plot
	"""
	def generateAttractorVsThreshold(self, thresholdMin, thresholdMax,
										thresholdStep, numICsPerThreshold):

		thetas = []
		attractors = []

		for t in np.arange(thresholdMin, thresholdMax, thresholdStep):
			thetas.append(t)

			attractors.append(self.countAttractors(numICs=numICsPerThreshold, threshold=t))	
			self.reset(resetEdges=False)

	    #Now, plot this using matplotlib
		plt.clf()
		plt.plot(thetas, attractors)
		plt.xlabel('Threshold')
		plt.ylabel('Attractor Count')	
		plt.grid(True)
		plt.savefig("")


	"""
	Create attractor count and GED vs threshold plot
	"""
	def generateDiversityVsThreshold(self, thresholdMin, thresholdMax,
									 thresholdStep, numICsPerThreshold):

		thetas = []
		attractors = []
		expressionDiversities = []

		for t in np.arange(thresholdMin, thresholdMax, thresholdStep):

			thetas.append(t)
			attractors.append(self.countAttractors(numICs=numICsPerThreshold, threshold=t))	
			expressionDiversities.append(self.expressionDiversity(threshold=t, numICs=numICsPerThreshold))	
			self.reset(resetEdges=False)

		#Normalize values
		maxA = max(attractors)
		maxE = max(expressionDiversities)

		attractorsPlot = [float(x) / maxA for x in attractors]
		expDivPlot = [float(x) / maxE for x in expressionDiversities]

	    #Now, plot this using matplotlib
		plt.clf()
		plt.plot(thetas, attractorsPlot,'g--^', label='Attractor Count')
		plt.plot(thetas, expDivPlot,'b-o', label='Expression Diversity')
		plt.legend( loc='upper left', numpoints = 1 )
		plt.xlabel('Threshold')
		plt.ylabel('Diversity Parameter')	
		plt.title('Diversity Parameter with varying threshold for 25 gene network')
		plt.grid(True)
		plt.savefig("AttractorGED.jpg")


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




		