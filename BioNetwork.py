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


	def get_degrees(self):
		degrees = []
		for i in xrange(0,self.n,1):
			degrees.append( sum(self.g[i]) )
		return degrees


	# Display network
	def __str__(self):
		# Matrix convertion to NetworkX graph:
		gx = networkx.to_networkx_graph(data=self.g)
		# Draw the graph with Matplotlib
		fig1 = pyplot.figure()
		networkx.draw(gx, with_labels=False, node_size=50, linewidths=0, width=0.5, alpha=1, cmap=pyplot.cm.rainbow, node_color=self.get_degrees())
		pyplot.savefig("NetworkX_plot1.png")
		# Draw the graph with Matplotlib
		fig2 = pyplot.figure()
		networkx.draw_networkx(gx, with_labels=False, node_size=50, linewidths=0, width=0.5, alpha=1, cmap=pyplot.cm.rainbow, node_color=self.get_degrees())
		pyplot.savefig("NetworkX_plot2-networkx.png")
		# Draw the graph with a circular layout.
		fig3 = pyplot.figure()
		networkx.draw_circular(gx, with_labels=False, node_size=50, linewidths=0, width=0.5, alpha=1, cmap=pyplot.cm.rainbow, node_color=self.get_degrees())
		pyplot.savefig("NetworkX_plot3-circular.png")
		# Draw the graph with a random layout.
		fig4 = pyplot.figure()
		networkx.draw_random(gx, with_labels=False, node_size=50, linewidths=0, width=0.5, alpha=1, cmap=pyplot.cm.rainbow, node_color=self.get_degrees())
		pyplot.savefig("NetworkX_plot4-random.png")
		# Draw the graph with a spectral layout.
		fig5 = pyplot.figure()
		networkx.draw_spectral(gx, with_labels=False, node_size=50, linewidths=0, width=0.5, alpha=1, cmap=pyplot.cm.rainbow, node_color=self.get_degrees())
		pyplot.savefig("NetworkX_plot5-spectral.png")
		# Draw the graph with a spring layout.
		fig6 = pyplot.figure()
		networkx.draw_spring(gx, with_labels=False, node_size=50, linewidths=0, width=0.5, alpha=1, cmap=pyplot.cm.rainbow, node_color=self.get_degrees())
		pyplot.savefig("NetworkX_plot6-spring.png")
		# Draw the graph with a shell layout.
		fig7 = pyplot.figure()
		networkx.draw_shell(gx, with_labels=False, node_size=50, linewidths=0, width=0.5, alpha=1, cmap=pyplot.cm.rainbow, node_color=self.get_degrees())
		pyplot.savefig("NetworkX_plot7-shell.png")
		# pyplot.draw()
		# pyplot.show()
		return "\nNetwork display saved.\n"





#============================================================================================================================
#						Main Script
#============================================================================================================================
def main():
	print "\n-----------------------------------------------------------------\n"

	# New network:
	n = Network(n=30)
	print n.g
	print n.get_degrees()
	print n
	

	print "\nExcecution successful."
	print "-----------------------------------------------------------------\n"

main()