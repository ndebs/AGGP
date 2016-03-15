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
		self.costSmallWorld=0
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
		# Test if the figure environment has already been opened:
		try:
			fig
		except: # No figure environment opened yet
			fig = pyplot.figure() # Opens a figure environment
		else: # Figure environment already yet
			pyplot.clf()
		# Draw the graph with Matplotlib
		networkx.draw(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot1.png")
		pyplot.clf()
		# Draw the graph with Matplotlib
		networkx.draw_networkx(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot2-networkx.png")
		pyplot.clf()
		# Draw the graph with a circular layout.
		networkx.draw_circular(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot3-circular.png")
		pyplot.clf()
		# Draw the graph with a random layout.
		networkx.draw_random(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot4-random.png")
		pyplot.clf()
		# Draw the graph with a spectral layout.
		networkx.draw_spectral(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot5-spectral.png")
		pyplot.clf()
		# Draw the graph with a spring layout.
		networkx.draw_spring(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot6-spring.png")
		pyplot.clf()
		# Draw the graph with a shell layout.
		networkx.draw_shell(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig("NetworkX_plot7-shell.png")
		pyplot.clf()
		# pyplot.draw()
		# pyplot.show()	
		# pyplot.close(fig)
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
			for j in xrange(i,self.n,1):
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
		self.dist = self.dist.astype(dtype=int)



	# Return the repartition vector of a noral path distribution
	def refPathRep(self,length):
		ref_rep = []
		sd = 1
		mean = math.log(math.log(self.n))
		for x in xrange(0,length,1):
			pnorm = math.exp(-(x-mean)**2/(2*sd**2))/(sd*math.sqrt(2*math.pi))
			ref_rep.append( pnorm*sum(range(0,self.n,1)) )
		# Complete the repartition until 0
		while ( pnorm*sum(range(0,self.n,1)) > 0.5 ):
			x += 1
			pnorm = math.exp(-(x-mean)**2/(2*sd**2))/(sd*math.sqrt(2*math.pi))
			ref_rep.append( pnorm*sum(range(0,self.n,1)) )
		return ref_rep



	# Return the network cost due to the small world constraint
	def smallWorldCost(self, plot=False):
		self.costSmallWorld = 0
		# Shortest paths between each pair of vertices: self.dist
		self.pairedShortestPaths()
		# Shortest paths distribution of the network:
		self_rep = self.matrixToDistribution(self.dist)
		# Creation of the reference repartition
		ref_rep = self.refPathRep(length=len(self_rep))
		# If self_rep is shorter than ref_rep
		while (len(ref_rep)>len(self_rep)):
			self_rep.append(0)
		# PRINTS
		"""
		print "Dist:\n",self.dist
		print "Shortest path distribution:\t",self_rep
		print "Shortest path normal distribution:\t",ref_rep
		"""
		# Computation of the Sum Square
		for i in xrange(0,len(self_rep),1):
			self.costSmallWorld += (self_rep[i] - ref_rep[i])**2
		if ( plot==True ):
			# Test if the figure environment has already been opened:
			try:
				fig
			except: # No figure environment opened yet
				fig = pyplot.figure() # Opens a figure environment
			else: # Figure environment already yet
				pyplot.clf()
			pyplot.plot(range(0,len(self_rep),1), self_rep, label='Network shortest path distribution', linestyle='--', marker='o', linewidth=1, markersize=5, color='red')
			pyplot.plot(range(0,len(self_rep),1), ref_rep, label='Normal shortest path distribution', linestyle='-', marker='.', linewidth=1, markersize=10, color='blue')
			# Plot Parameters
			pyplot.xlabel("Shortest path length")
			pyplot.ylabel("Number")
			pyplot.legend(fontsize=10) #adds a legend
			pyplot.show()
			pyplot.clf()
			pyplot.close(fig)





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
		# return the 3 meands and variances of the costs, output data are arranged in alphabetical order
		meanCostClique=0
		varCostClique=0
		meanCostSmallWorld=0
		varCostSmallWorld=0
		meanCostPowerLaw=0
		varCostPowerLaw=0
		for i,e in enumerate(self.graphs):
			#Calculation for cliques
			meanCostClique+=e.costClique
			varCostClique+=math.pow(e.costClique,2)
			#Calculation for small world
			meanCostSmallWorld+=e.costSmallWorld
			varCostSmallWorld+=math.pow(e.costSmallWorld,2)
			#Calculation for power law
			# COMING SOON
		meanCostClique=meanCostClique/self.m
		varCostClique=(varCostClique/self.m)-math.pow(meanCostClique,2)
		meanCostSmallWorld=meanCostSmallWorld/self.m
		varCostSmallWorld=(varCostSmallWorld/self.m)-math.pow(meanCostSmallWorld,2)
		return [meanCostClique, math.sqrt(varCostClique), meanCostSmallWorld, math.sqrt(varCostSmallWorld)]

	def overallCost(self):
		sdCost=self.sdPopCost()
		for i,e in enumerate(self.graphs):
			#Normalization of the 3 costs
			e.cost=(e.costClique-sdCost[0])/sdCost[1]+(e.costSmallWorld-sdCost[2])/sdCost[3] # + others costs

	def averagePopCost(self):
		# return the average cost of networks in population
		averageCost=0
		for i,e in enumerate(self.graphs):
			averageCost=averageCost+e.cost
		#print averageCost
		averageCost=averageCost/float(self.m)
		print "Average cost in population:"
		print averageCost
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
	n.smallWorldCost(plot=True)
	print "Small World Cost =\t",n.costSmallWorld
	
	m=5
	nodes=10
	print "\nNetwork population: %d graphs with %d nodes " %(m,nodes)
	P=Population(m,nodes)
	for i in xrange(m):
		P.graphs[i].cliqueCost()
		P.graphs[i].smallWorldCost()
	
	print "Standard deviance  in population: %f, %f" % (P.sdPopCost()[0], P.sdPopCost()[1])
	P.overallCost()
	P.averagePopCost()


	print "\nExcecution successful."
	print "-----------------------------------------------------------------\n"

main()