#############################################################
#							AGGP							#
#				Biological Network generator				#
#	Noelie		Marianne		Nathalie		Vincent		#
#############################################################
import matplotlib.pyplot as pyplot
import numpy as np
import math 
import networkx as networkx



#============================================================================================================================
#					Class NETWORK
#============================================================================================================================
class Network(object):

	# Constructor
	def __init__(self,n=1):
		self.n = n
		self.g = np.zeros(shape=(self.n,self.n), dtype=np.int8)
		self.cost=0 # sum of the 3 costs
		self.costClique=0
		for i in xrange(0,self.n,1):
			for j in xrange(i+1,self.n,1):
				self.g[i,j] = np.random.random_integers(low=0, high=1, size=None)
				self.g[j,i] = self.g[i,j] # Symmetrical matrix



	# Display network
	def __str__(self):
		# Matrix convertion to NetworkX graph:
		gx = networkx.to_networkx_graph(data=self.g)
		# Parameters
		deg = self.get_degrees()
		size_deg = [(10+(90*(i-min(deg)))/float(max(deg))) for i in deg] # size_deg = [(100*i/float(max(deg))) for i in deg]
		numlabels = False
		edgecolor = 'grey'
		nodecmap = pyplot.cm.rainbow
		# Draw the graph with Matplotlib
		fig1 = pyplot.figure()
		networkx.draw(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot1.png")
		# Draw the graph with Matplotlib
		fig2 = pyplot.figure()
		networkx.draw_networkx(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot2-networkx.png")
		# Draw the graph with a circular layout.
		fig3 = pyplot.figure()
		networkx.draw_circular(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot3-circular.png")
		# Draw the graph with a random layout.
		fig4 = pyplot.figure()
		networkx.draw_random(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot4-random.png")
		# Draw the graph with a spectral layout.
		fig5 = pyplot.figure()
		networkx.draw_spectral(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot5-spectral.png")
		# Draw the graph with a spring layout.
		fig6 = pyplot.figure()
		networkx.draw_spring(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot6-spring.png")
		# Draw the graph with a shell layout.
		fig7 = pyplot.figure()
		networkx.draw_shell(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot7-shell.png")
		# pyplot.draw()
		# pyplot.show()
		return "\nNetwork display saved.\n"



	def get_degrees(self):
		degrees = []
		for i in xrange(0,self.n,1):
			degrees.append( sum(self.g[i]) )
		return degrees



	def cliqueCost(self,coeffA=1,coeffB=1):
		l=self.g.shape[0]
		k=np.zeros(l)
		c=np.zeros(l)
		i=0
		while i<l:
			j=0
			while j<l:
				if self.g[i,j]==1:
					k[i]+=1
					m=j
					while m<l:
						if self.g[i,m]==1 and self.g[j,m]==1:
							c[i]+=1
						m+=1
				j+=1
			if k[i]>1: # avoid to divide by 0
				c[i]=(2*c[i])/(k[i]*(k[i]-1))
			i+=1
		a,b=np.polyfit(k,c,1) # polynomial regression, degree one
		self.costClique=coeffA*abs(a)+coeffB*(b-np.mean(c))
	


	# Return the value repartition of a matrix
	def matrixToDistribution(self,mat):
		# Headcount initialized to 0
		rep = [0]*(np.amax(a=mat)+1)
		# Compute the distribution
		for i in xrange(0,self.n,1):
			for j in xrange(i+1,self.n,1):
				rep[int(mat[i,j])] += 1
		return rep



	# Return the matrix of minimum distances between each pair of nodes
	def pairedShortestPaths(self):
		self.dist = np.empty(shape=(self.n,self.n))
		self.dist[:] = float('Inf')
		# Set path to 0 for each node
		for i in xrange(0,self.n,1):
			self.dist [i,i] = 0
		# Set path to 1 for each edge
		for i in xrange(0,self.n,1):
			for j in xrange(i+1,self.n,1):
				if (self.g[i,j] == 1):
					self.dist[i,j] = 1
					self.dist[j,i] = 1
		# Compute the shortest path for each pair of nodes
		for k in xrange(0,self.n,1):
			for i in xrange(0,self.n,1):
				for j in xrange(0,self.n,1):
					if (self.dist[i,j] > self.dist[i,k]+self.dist[k,j]):
						self.dist[i,j] = self.dist[i,k]+self.dist[k,j]
		return self.dist.astype(dtype=int)



#============================================================================================================================
#					Class POPULATION
#============================================================================================================================

class Population(object):

	def __init__(self,m,n):
		self.m=m # number of networks in population
		self.n=n # node number in one network
		self.graphs=range(m)
		for i in xrange(m):
			self.graphs[i]=Network(n=n)


	def sdPopCost(self):
		# Have to add the 2 othets costs
		meanCostClique=0
		varCostClique=0
		for i,e in enumerate(self.graphs):
			meanCostClique+=e.costClique
			varCostClique+=math.pow(e.costClique,2)
		meanCostClique=meanCostClique/self.m
		varCostClique=(varCostClique/self.m)-math.pow(meanCostClique,2)
		return math.sqrt(varCostClique)

	def overallCost(self):
		for i,e in enumerate(self.graphs):
			e.cost=e.costClique # + others costs

	def averagePopCost(self):
		# return the average cost of networks in population
		averageCost=0
		for i,e in enumerate(self.graphs):
			averageCost+=e.cost
		averageCost=averageCost/self.m
		return averageCost



#============================================================================================================================
#						Main Script
#============================================================================================================================

def main():
	print "\n-----------------------------------------------------------------\n"

	# New network:
	n = Network(n=15)
	print n.g
	print n.get_degrees()
	print n
	n.cliqueCost(1,1)
	print "Dist:\n",n.pairedShortestPaths()
	print n.matrixToDistribution(n.dist)
	
	m=2
	nodes=10
	print "\n \n Network population: %d graphs with %d nodes " %(nodes,m)
	P=Population(m,nodes)
	for i in xrange(m):
		P.graphs[i].cliqueCost()
		print P.graphs[i].costClique
	print "Standard deviance  in population: %f" % P.sdPopCost()
	P.overallCost()
	print "Average cost in population: %f " % P.averagePopCost()
	print "\nExcecution successful."
	print "-----------------------------------------------------------------\n"

main()