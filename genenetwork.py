import scipy as sp
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import math
import time
import sys

"""""""""""""""""""""""""""""""""""""""""""""""
A geneNetwork class for gene network analysis.

James Komianos
Last Edited: 11/17/13
"""""""""""""""""""""""""""""""""""""""""""""""

class geneNetwork:

	G = nx.Graph()

	"""
	Constructor
	nodes= a list tuples containing node name and value
	edgeTuples= a list of tuples containing
		the edges to add, with weights
	"""
	def __init__(self, nodeTuples, edgeTuples):

		for nodeTuple in nodeTuples:
			self.G.add_node(nodeTuple[0], value=nodeTuple[1])

		self.G.add_weighted_edges_from(edgeTuples)

	"""
	Update graph using thresholding
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
				sumN += value * weight - threshold

			#Update node based on sign of sum	
			if (sumN >= 0): self.G.node[n]['value'] = 1 
			else: self.G.node[n]['value'] = 0

	"""
	For debugging purposes
	"""
	def printNetwork(self, outputFile):

		plt.clf()
		pos=nx.spring_layout(self.G)

		val_map = {1: 0.2564, 0: 1.0}
		values = [val_map.get(self.G.node[n]['value'], 0.1) for n in self.G]

		nx.draw(self.G, pos, cmap=plt.get_cmap('prism'), node_color=values)
		plt.savefig(outputFile)




		