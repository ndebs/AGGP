#############################################################
#							AGGP							#
#				Biological Network generator				#
#	Noelie		Marianne		Nathalie		Vincent		#
#############################################################
import matplotlib.pyplot as pyplot
import numpy as np
import math 
import networkx as networkx
import random
import scipy.stats
import copy as copy

#============================================================================================================================
#					Class NETWORK
#============================================================================================================================
class Network(object):

	# Constructor
	def __init__(self,n=1):
		self.n = n
		self.g = np.zeros(shape=(self.n,self.n), dtype=np.int8)
		self.cost=0 # sum of the 3 costs
		self.costRelative=0 # sum of the 3 costs
		self.costClique=0
		self.costSmallWorld=0
		self.costDegree=0
		for i in xrange(0,self.n,1):
			for j in xrange(i+1,self.n,1):
				self.g[i,j] = np.random.random_integers(low=0, high=1, size=None)
				self.g[j,i] = self.g[i,j] # Symmetrical matrix
			# Prevent the construction of a graph with a non-connected node
			if (sum(self.g[i])==0):
				# Additon of a random vertex
				r = np.random.random_integers(low=0, high=self.n-1, size=None)
				while (i==r):
					r = np.random.random_integers(low=0, high=self.n-1, size=None)
				self.g[i,r] = 1
				self.g[r,i] = self.g[i,r] # Symmetrical matrix



	# Display network
	def __str__(self,ID=''):
		# Matrix convertion to NetworkX graph:
		gx = networkx.to_networkx_graph(data=self.g)
		# Parameters
		deg = self.get_degrees()
		size_deg = [(10+(90*(i-min(deg)))/float(max(deg))) for i in deg] # size_deg = [(100*i/float(max(deg))) for i in deg]
		numlabels = False
		edgecolor = 'grey'
		nodecmap = pyplot.cm.rainbow
		# # Draw the graph with Matplotlib
		# networkx.draw(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		# pyplot.savefig(ID+"NetworkX_plot1.png")
		# pyplot.clf()
		# # Draw the graph with Matplotlib
		# networkx.draw_networkx(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		# pyplot.savefig(ID+"NetworkX_plot2-networkx.png")
		# pyplot.clf()
		# Draw the graph with a circular layout.
		networkx.draw_circular(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		pyplot.savefig(ID+"NetworkX_plot3-circular.png")
		pyplot.clf()
		# # Draw the graph with a random layout.
		# networkx.draw_random(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		# pyplot.savefig(ID+"NetworkX_plot4-random.png")
		# pyplot.clf()
		# # Draw the graph with a spectral layout.
		# networkx.draw_spectral(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		# pyplot.savefig(ID+"NetworkX_plot5-spectral.png")
		# pyplot.clf()
		# # Draw the graph with a spring layout.
		# networkx.draw_spring(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		# pyplot.savefig(ID+"NetworkX_plot6-spring.png")
		# pyplot.clf()
		# # Draw the graph with a shell layout.
		# networkx.draw_shell(gx, with_labels=numlabels, node_size=size_deg, linewidths=0, width=0.5, alpha=1, cmap=nodecmap, node_color=deg, edge_color=edgecolor)
		# pyplot.savefig(ID+"NetworkX_plot7-shell.png")
		# pyplot.clf()
		# pyplot.draw()
		# pyplot.show()
		return "\nNetwork display saved.\n"

	def fileCytoscape(self):
		f=open('graphCytoscape.txt', 'w')
		for i in xrange(0,self.n,1):
			for j in xrange(i+1,self.n,1):
				if self.g[i,j]==1:
					f.write('%d %d\n' % (i,j))
		f.close()




	def get_degrees(self):
		degrees = []
		for i in xrange(0,self.n,1):
			degrees.append( sum(self.g[i]) )
		return degrees


	def mutation(self,mut_rate):
		for i in xrange(0,self.n,1):
			for j in xrange(i+1,self.n,1):
				rand_freq = random.uniform(0,1)
				if rand_freq < mut_rate:
					if self.g[i,j] == 1:
						self.g[i,j] = 0
						self.g[j,i] = self.g[i,j]
					if self.g[i,j] == 0:
						self.g[i,j] = 1
						self.g[j,i] =  self.g[i,j]
			# Prevent the mutation of a graph with a non-connected node
			if (sum(self.g[i])==0):
				# Additon of a random vertex
				r = np.random.random_integers(low=0, high=self.n-1, size=None)
				while (i==r):
					r = np.random.random_integers(low=0, high=self.n-1, size=None)
				self.g[i,r] = 1
				self.g[r,i] = self.g[i,r] # Symmetrical matrix



############################################################################################################################################
#                              CONSTRAINT CLIQUE
############################################################################################################################################

	def cliqueCost(self):
		self.costClique=0
		l=self.g.shape[0]
		k=np.zeros(l)
		c=np.zeros(l)
		log_k=np.zeros(l)
		log_c=np.zeros(l)
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
			if c[i]<=0:
				c[i]=0.0000000001
			log_k[i]=math.log(k[i],10)
			log_c[i]=math.log(c[i],10)
			i+=1

		a,b=np.polyfit(log_k,log_c,1) # polynomial regression, degree one, older version
		for i, e in enumerate(log_k):
			self.costClique+=(a*e-e)**2 + 0.05*b*0 # b is minor in the calculation of cost
		


############################################################################################################################################
#                              CONSTRAINT SMALL WORLD
############################################################################################################################################

	# Return the value repartition of a matrix
	def matrixToDistribution(self,mat):
		# Headcount initialized to 0
		rep = [0]*(np.amax(a=mat)+1)
		# Compute the distribution
		for i in xrange(0,self.n,1):
			for j in xrange(i,self.n,1):
				# try:
				# 	rep[int(mat[i,j])]
				# except: 
				# 	print self.g
				# 	print self
				# 	print mat
				# 	print "index=",int(mat[i,j]),"\t len(rep)=",len(rep),"\t rep=",rep
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
		# Shortest paths between each pair of vertices: self.dist
		self.pairedShortestPaths()
		averageShortestPaths = np.average(a=self.dist)
		# Reference mean:
		meanRef = math.log(math.log(self.n))
		# Cost:
		self.costSmallWorld = math.fabs(averageShortestPaths - meanRef)

		# # # Shortest paths distribution of the network:
		# # self_rep = self.matrixToDistribution(self.dist)
		# # # Creation of the reference repartition
		# # ref_rep = self.refPathRep(length=len(self_rep))
		# # # If self_rep is shorter than ref_rep
		# # while (len(ref_rep)>len(self_rep)):
		# # 	self_rep.append(0)
		# # # PRINTS
		# # """
		# # print "Dist:\n",self.dist
		# # print "Shortest path distribution:\t",self_rep
		# # print "Shortest path normal distribution:\t",ref_rep
		# # """
		# # # Computation of the Sum Square
		# # for i in xrange(0,len(self_rep),1):
		# # 	self.costSmallWorld += (self_rep[i] - ref_rep[i])**2
		if ( plot==True ):
			print "Shortest path between each node:\n",self.dist
			print "Average shortest path:",averageShortestPaths
			print "log(log(n)) =",meanRef
			print "Small World Cost =",self.costSmallWorld
			# # Test if the figure environment has already been opened:
			# # try:
			# # 	fig
			# # except: # No figure environment opened yet
			# # 	fig = pyplot.figure() # Opens a figure environment
			# # else: # Figure environment already yet
			# # 	pyplot.clf()
			# pyplot.plot(range(0,len(self_rep),1), self_rep, label='Network shortest path distribution', linestyle='--', marker='o', linewidth=1, markersize=5, color='red')
			# pyplot.plot(range(0,len(self_rep),1), ref_rep, label='Normal shortest path distribution', linestyle='-', marker='.', linewidth=1, markersize=10, color='blue')
			# # Plot Parameters
			# pyplot.xlabel("Shortest path length")
			# pyplot.ylabel("Number")
			# pyplot.legend(fontsize=10) #adds a legend
			# pyplot.show()
			# pyplot.clf()


############################################################################################################################################
#                              CONSTRAINT POWER DEGREE
############################################################################################################################################

	def P(self,k,gamma):
		return k**(-gamma)

	
	def P_obs(self):
		x_list=range(1,20)
		l = self.get_degrees() #list of all degrees (for each node)
		freq_list= []
		for i in x_list:
			ni= l.count(i) #ni = number of nodes of degree i
			freq_list.append(float(ni)/float(self.n))
		return (x_list,freq_list)



	def degreeCost(self,gamma):
		self.costDegree = 0 
		degrees,freq_list = self.P_obs()
		for i in xrange(len(freq_list)-1): 
			#compute cost (=SCE)
			fk= freq_list[i]
			k= degrees[i]
			self.costDegree += (fk - self.P(k,gamma))**2

	def plot_freq_degree(self,gamma_opti):
		y_theo= []
		y_obs= []
		k_list=[]
		x=xrange(1,20)
		for k in self.get_degrees(): #k = degree of each vertex
			k_list.append(k)
			fk=  float(k)/float(self.n)
			y_obs.append(fk)
		for i in x:
			y_theo.append(self.P(i,gamma_opti)) 
		#print y_obs
		pyplot.plot(x,y_theo,label='Theoric degree distribution', linestyle='--', marker='o', linewidth=1, markersize=5, color='red')
		#pyplot.plot(k_list,y_obs, label='Observed distribution', marker='+', linewidth=1, markersize=10, color='blue')
		pyplot.hist(k_list,normed=1,label='Observed distribution')
		pyplot.xlabel("degree k")
		pyplot.ylabel("Frequency")
		pyplot.show()
'''
	def d_P(self,k,gamma):
		return -gamma*k**(-gamma-1)

	def d_degreeCost(self,gamma):
		d = 0 
		degrees,freq_list = self.P_obs()
		for i in xrange(len(freq_list)-1): 
			fk= freq_list[i]
			k= degrees[i]
			d += (fk - self.P(k,gamma))*self.d_P(k,gamma)
		return -d
		

	def d2_degreeCost(self,gamma):
		d2 = 0 
		degrees,freq_list = self.P_obs()
		for i in xrange(len(freq_list)-1): 
			k= degrees[i]
			d2 += self.d_P(k,gamma)**2 
		return d2
		

	def algo_LM(self,maxit):
		#initialisation of at 2
		gamma = 2
		#definition of a parameter lambda for LM algo
		lamb = 0.001
		self.degreeCost(gamma)
		cost = self.costDegree
		# defintion of the vector d_degreeCost's norm
		norme_cost = [abs(self.d_degreeCost(gamma))]
		#iteration until maximum iteration maxit
		i=0
		while (i<maxit):
			Grad_gamma = self.d_degreeCost(gamma)
			H_LM = self.d2_degreeCost(gamma)*(1+lamb)
			d_LM = -Grad_gamma/H_LM
			self.degreeCost(gamma+d_LM)
			new_cost = self.costDegree
			#print "cost","new_cost", cost, new_cost
			if new_cost<cost:
			#change gamma
				cost = new_cost
				gamma = gamma+d_LM
				lamb=lamb/10
				norme_cost = [abs(self.d_degreeCost(gamma))]
			else:
				#don't change gamma
				lamb=lamb*10
			i += 1
		#return gamma optimal
		return gamma
'''




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
			meanCostPowerLaw+=e.costDegree
			varCostPowerLaw+=math.pow(e.costDegree,2)

		meanCostClique=meanCostClique/float(self.m)
		varCostClique=(varCostClique/float(self.m))-math.pow(meanCostClique,2)
		meanCostSmallWorld=meanCostSmallWorld/float(self.m)
		varCostSmallWorld=(varCostSmallWorld/float(self.m))-math.pow(meanCostSmallWorld,2)
		meanCostPowerLaw=meanCostPowerLaw/float(self.m)
		varCostPowerLaw=(varCostPowerLaw/float(self.m))-math.pow(meanCostPowerLaw,2)
		return [meanCostClique, math.sqrt(math.fabs(varCostClique)), meanCostSmallWorld, math.sqrt(math.fabs(varCostSmallWorld)), meanCostPowerLaw, math.sqrt(math.fabs(varCostPowerLaw))]

	def overallCost(self):#,coeffCli=1,coeffSwl=20,coeffDeg=5):
		sdCost=self.sdPopCost()
		costPerGraph=[]
		costPerGraphClique=[]
		costPerGraphSwallWorld=[]
		costPerGraphDegree=[]
		costPerGraphRelative=[]
		# print "SD",sdCost
		# Save costs of the population in vectors
		for i,e in enumerate(self.graphs):
			costPerGraphClique.append(e.costClique)
			costPerGraphSwallWorld.append(e.costSmallWorld)
			costPerGraphDegree.append(e.costDegree)
			print "Clique = ",e.costClique,"\tSW =",e.costSmallWorld, "\tDeg =",e.costDegree 
		maxCosts = [max(costPerGraphClique),max(costPerGraphSwallWorld),max(costPerGraphDegree)]
		print "maxCosts = ",maxCosts
		"""for i,e in enumerate(self.graphs): # Normalization of the 3 costs and sum of them
			epsilon = 10**(-10)
			# e.cost=(e.costClique-sdCost[0])/float(sdCost[1]+epsilon)+(e.costSmallWorld-sdCost[2])/float(sdCost[3]+epsilon) + (e.costDegree-sdCost[4])/float(sdCost[5]+epsilon)
			e.cost= coeffCli*e.costClique + coeffSwl*e.costSmallWorld + coeffDeg*e.costDegree
			print "Clique = ",coeffCli*e.costClique/e.cost,"\tSW =", coeffSwl*e.costSmallWorld/e.cost, "\tDeg =", coeffDeg*e.costDegree/e.cost
			costPerGraph.append(e.cost)"""
		for i,e in enumerate(self.graphs):
			# Absolute costs (used for plotting)
			e.cost = costPerGraphClique[i]+costPerGraphSwallWorld[i]+costPerGraphDegree[i]
			# print "Graph cost", e.cost
			costPerGraph.append(e.cost)
			# Relative costs
			e.costRelative = costPerGraphClique[i]/float(3*maxCosts[0])+costPerGraphSwallWorld[i]/float(3*maxCosts[1])+costPerGraphDegree[i]/float(3*maxCosts[2])
			costPerGraphRelative.append(e.costRelative)
		print "Absolute cost per graph = ",costPerGraph
		print "Relative cost per graph = ",costPerGraphRelative
		return costPerGraph,costPerGraphRelative

	def averagePopCost(self):
		# return the average cost of networks in population
		averageCost=0
		for i,e in enumerate(self.graphs):
			averageCost+=e.cost
		averageCost=averageCost/float(self.m)
		print "Average cost in population:"
		print averageCost
		return averageCost


	def selection(self, costPerGraph, c=0.5):
		fitness = [1-i for i in costPerGraph]
		rank = scipy.stats.rankdata(fitness)
		print "Fitness:",fitness,"\nRank:",rank
		# ELITISM
		# Initiate a new pop with the best graph
		print "Argmax =",np.argmax(fitness)
		newPop = [ copy.deepcopy(self.graphs[np.argmax(fitness)]) ]
		# print "BEST:\tr=",rank[np.argmax(fitness)],"\tfit=",fitness[np.argmax(fitness)]
		# RANKING
		Wr = [] # Proba of reproducing (non-normalized)
		for i,e in enumerate(self.graphs):
			# Reproduction probability of each graph:
			Wr.append( self.m*(c-1)*(c**(self.m-rank[i]))/(c**(self.m)-1) )
			# print "r=",rank[i],"\tfit=",fitness[i],"\tProba = ",Wr
		Wr = [i/float(sum(Wr)) for i in Wr]
		# print "Wr = ",Wr
		nbRepro = np.random.multinomial(n=(self.m-1), pvals=Wr, size=None)
		print "Offspring: ",nbRepro
		for i,e in enumerate(self.graphs):
			for j in xrange(1,nbRepro[i]+1,1):
				newPop.append( copy.deepcopy(e) )
		# for g in newPop:
		# 	print g.cost,
		# print "\n"
		self.graphs = copy.deepcopy(newPop)


	def crossingOver(self, tx=0.05):
		p=np.random.binomial(self.m, tx) # p crossing overs have to be made
		#print p
		while p>0:
			i=random.randint(1,self.m-1)
			j=random.randint(1,self.m-1)
			while i==j:
				j=random.randint(1,self.m-1)
			ncol=random.randint(0,self.n-1) # row/colum to change 
			#print ncol
			tempI=self.graphs[i].g[0:self.n,ncol] # intermediar copy is mandatory
			temp_J=self.graphs[j].g[0:self.n,ncol]
			self.graphs[i].g[ncol,0:self.n]=temp_J
			self.graphs[j].g[ncol,0:self.n]=tempI
			self.graphs[i].g[0:self.n,ncol]=self.graphs[i].g[ncol,0:self.n]
			self.graphs[j].g[0:self.n,ncol]=self.graphs[j].g[ncol,0:self.n]
			p-=1



	def updatePop(self,generation,gamma,c,mut_rate,cro_rate):
		# Variable to keep a trace of:
		minGenerationCost = []
		bestGenerationCost = []
		aveGenerationCost = []
		# Population evolution
		for generationNumber in xrange(0,generation,1):
			# Initial costs
			# print "Initial cost:",
			for i,G in enumerate(self.graphs):
				G.cliqueCost()
				G.smallWorldCost()
				G.degreeCost(gamma=gamma)
				# print G.cost,
			print "\n------------------------------------\nNEW GENERATION:\t"+str(generationNumber)+"\n-------------------------------------\n"
			# Save population costs
			popAboluteCost,popRelativeCost = self.overallCost() # ABOLUTE,RELATIVE
			print "Min abs = ",min(popAboluteCost)
			print "Best = ",self.graphs[0].cost
			minGenerationCost.append(min(popAboluteCost))
			bestGenerationCost.append(self.graphs[0].cost)
			aveGenerationCost.append(self.averagePopCost())
			# Selection
			self.selection(costPerGraph=popRelativeCost, c=c)
			# Mutation of each graph (except the best one)
			for i in xrange(1,self.m,1):
				self.graphs[i].mutation(mut_rate=mut_rate)
			# Crossing over
			self.crossingOver(tx=cro_rate)
			self.graphs[0].__str__(ID="0-"+str(generationNumber)+"-")
			self.graphs[1].__str__(ID="1-"+str(generationNumber)+"-")
		self.graphs[0].fileCytoscape()
		# Costs plot:
		pyplot.plot(range(0,generation,1), minGenerationCost, label='Minimum', linestyle='-', marker='.', linewidth=1, markersize=5, color='red')
		pyplot.plot(range(0,generation,1), bestGenerationCost, label='Best conserved', linestyle='--', marker='.', linewidth=1, markersize=5, color='green')
		pyplot.plot(range(0,generation,1), aveGenerationCost, label='Average', linestyle='--', marker='.', linewidth=1, markersize=5, color='blue')
		# Plot Parameters
		pyplot.xlabel("Shortest path length")
		pyplot.ylabel("Number")
		pyplot.legend(fontsize=10) #adds a legend
		pyplot.show()





#============================================================================================================================
#						Main Script
#============================================================================================================================

def main():
	print "\n-----------------------------------------------------------------\n"

	# New network:
	
	n = Network(n=15)
	print n.g
	print n.get_degrees()
	# print n

	# CLIQUE COST
	n.cliqueCost()
	print "Clique cost ", n.costClique
	# SMALL WORD
	n.smallWorldCost(plot=True)
	# print "Small World Cost =\t",n.costSmallWorld
	# POWER DEGREE
	gamma=2.2
	n.degreeCost(gamma) 
	print "cost degree", n.costDegree
	# n.plot_freq_degree(gamma_opti)
	n.fileCytoscape()
	print n

	m=5
	nodes=10
	print "\nNetwork population: %d graphs with %d nodes " %(m,nodes)
	P=Population(m,nodes)
	
	for i,G in enumerate(P.graphs):
		G.cliqueCost()
		G.smallWorldCost()
		G.degreeCost(gamma)
	
	print "Standard deviance  in population: %f, %f" % (P.sdPopCost()[0], P.sdPopCost()[1])
	b=P.overallCost()
	a=P.averagePopCost()

	print "Test crossing over"
	P_cross=Population(2,4)
	print P_cross.graphs[0].g
	print P_cross.graphs[1].g
	P_cross.crossingOver()
	print P_cross.graphs[0].g
	print P_cross.graphs[1].g
	P.crossingOver()

	print "\nExcecution successful."
	print "-----------------------------------------------------------------\n"

#main()

P=Population(m=10,n=25)
P.updatePop(generation=10,gamma=2.2,c=0.5,mut_rate=0.15,cro_rate=0.05)