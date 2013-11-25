import scipy as sp
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import math
import time
import sys

"""""""""""""""""""""""""""""""""""""""""""""""
A geneNetwork class for gene network analysis.

James Komianos
Last Edited: 11/20/13
"""""""""""""""""""""""""""""""""""""""""""""""

class geneNetwork:

	G = nx.Graph()
	GPerm = None
	numGenes = 0

	"""
	Constructor
	nodes= a list tuples containing node name and value
	edgeTuples= a list of tuples containing
		the edges to add, with weights
	"""
	def __init__(self, preset="RANDOM", numNodes=100, 
				connectionsPerNode=5, nodeTuples=None, edgeTuples=None):

		if preset == "RANDOM":
			#Add nodes
			for i in xrange(0, numNodes):

				state = 1 if random.random() >= 0.5 else 0
				self.G.add_node(i, value=state)

			#Now add edges
			for i in xrange(0, numNodes):
				for j in xrange(0, connectionsPerNode):
					randomNode = random.randint(0,numNodes-1)
					randomWeight = 1 if random.random() >= 0.25 else -1
					toConnect = randomNode if randomNode != i else random.randint(0,numNodes-1)
					self.G.add_edge(i, randomNode, weight=randomWeight)

			self.numGenes = numNodes

		if preset == "MANUAL":
			for nodeTuple in nodeTuples:
				self.G.add_node(nodeTuple[0], value=nodeTuple[1])

			self.G.add_weighted_edges_from(edgeTuples)
			self.numGenes = len(nodeTuples)

	"""
	Update graphs using thresholding
	"""	
	def update(self, threshold):

		#Iterate over all nodes
		for n in self.G:
			sumN = 0
			value = self.G.node[n]['value']
			#Find all adjacents
			for neighbor in self.G.neighbors_iter(n):
				#Get edge weight
				weight = self.G[n][neighbor]['weight']
				sumN += value * weight

			#Update node based on sign of sum	
			print "Node=", n, "Value of change", sumN - threshold, "Current Value=",value
			if (sumN - threshold >= 0): self.G.node[n]['value'] = 1
			else: self.G.node[n]['value'] = 0

		#Update the permutation if it exists
		if self.GPerm != None:	
			for n in self.GPerm:
				sumN = 0
				value = self.GPerm.node[n]['value']
				#Find all adjacents
				for neighbor in self.GPerm.neighbors_iter(n):
					#Get edge weight
					weight = self.GPerm[n][neighbor]['weight']
					sumN += value * weight

				#Update node based on sign of sum	
				if (sumN - threshold>= 0): self.GPerm.node[n]['value'] = 1 
				else: self.GPerm.node[n]['value'] = 0

	"""
	Create permutation of current graph
	"""
	def createPermutation(self, diffThreshold=0.1):
		self.GPerm = self.G.copy()
		#Change percentage of nodes (keep edges the same)
		numToChange = int(diffThreshold * self.numGenes)

		for i in xrange(0, numToChange):
			nodeToChange = random.randint(0,self.numGenes - 1)
			currentVal = self.GPerm.node[nodeToChange]['value']
			self.GPerm.node[nodeToChange]['value'] = 1 if currentVal == 0 else 0

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
	def calculateSatHammingDist(self, numTimeSteps=20, threshold=0.5):	

		hammingDists = []
		for t in xrange(0, numTimeSteps):
			#Update, calculate hamming distance
			self.update(threshold)
			currentHammingDist = self.calculateHammingDist()
			hammingDists.append(currentHammingDist)

		return hammingDists	

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

		val_map = {1: 0.2564, 0: 1.0}
		values = [val_map.get(self.G.node[n]['value'], 0.1) for n in self.G]

		nx.draw(self.G, pos, cmap=plt.get_cmap('prism'), node_color=values)
		plt.savefig(outputFile)

	def printState(self):
		print "PRINTING CURRENT NODE DATA"
		for i in xrange(0, self.numGenes):
			print "Node=", i, "Value G=", self.G.node[i]['value'], "Value GPerm=", self.GPerm.node[i]['value']




		