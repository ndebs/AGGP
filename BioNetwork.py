#############################################################
#							AGGP							#
#				Biological Network generator				#
#	Noelie		Marianne		Nathalie		Vincent		#
#############################################################
import matplotlib.pyplot as pyplot
import numpy as numpy
import math as maths
import networkx as networkx



#============================================================================================================================
#					Class NETWORK
#============================================================================================================================
class Network(object):

	# Constructor
	def __init__(self,n=1,randomly=True):
		self.n = n
		self.g = numpy.zeros(shape=(self.n,self.n), dtype=numpy.int8)
		for i in xrange(0,self.n,1):
			for j in xrange(i+1,self.n,1):
				self.g[i,j] = numpy.random.random_integers(low=0, high=1, size=None)
				self.g[j,i] = self.g[i,j] # Symmetrical matrix



	# Display network
	def __str__(self):
		# Matrix convertion to NetworkX graph:
		gx = networkx.to_networkx_graph(data=self.g)
		# Draw
		networkx.draw(gx)
		pyplot.savefig("NetworkX_plot1.png")
		pyplot.show()
		# Draw
		networkx.draw_networkx(gx)
		pyplot.savefig("NetworkX_plot2-networkx.png")
		pyplot.show()
		# Draw circular:
		networkx.draw_circular(gx)
		pyplot.savefig("NetworkX_plot3-circular.png")
		pyplot.show()
		# Draw random ?
		networkx.draw_random(gx)
		pyplot.savefig("NetworkX_plot4-random.png")
		pyplot.show()
		# Draw random ?
		networkx.draw_spectral(gx)
		pyplot.savefig("NetworkX_plot5-spectral.png")
		# pyplot.draw()
		pyplot.show()
		# pyplot.draw()
		# pyplot.show()
		return "\n>>>\tNetwork displayed.\n"





#============================================================================================================================
#						Main Script
#============================================================================================================================
def main():
	print "\n-----------------------------------------------------------------\n"

	# New network:
	n = Network(n=15)
	print n.g
	print n
	

	print "\nExcecution successful."
	print "-----------------------------------------------------------------\n"

main()