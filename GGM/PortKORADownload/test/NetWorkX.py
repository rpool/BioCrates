#! /usr/bin/env python
# (F)utilities to explore the network of Escherichia coli
# as given on RegulonDB


##### FETCH THE NETWORK ON THE INTERNET ################################


import networkx as nx
from urllib2 import urlopen

# Fetch E. coli's network on Regulon DB as a .txt file

url = "http://regulondb.ccg.unam.mx/data/network_tf_gene.txt"
regulonFile = urlopen(url)

# Create an (empty) directed graph

eColiNetwork = nx.DiGraph(name = 'E. Coli Network')

# Read the file and fill our network

for line in regulonFile :

    # ignore comments and empty lines

    if not line.startswith(('\n','\t','#')):

        # A line is of the form "gene1  gene2   + "

        g1,g2, sign = line.split('\t')[:3]

        # standardize gene names ( all characters to lower case )

        g1,g2 = g1.lower(),g2.lower()

        # add the edge to the network

        eColiNetwork.add_edge(g1,g2)
        eColiNetwork[g1][g2]['sign'] = sign


########### SAVE AND LOAD ##############################################


import pickle

# Save the graph as EcoliNetwork.gn

with open('EcoliNetwork.gn', 'wb') as output:
    pickle.dump(eColiNetwork, output, pickle.HIGHEST_PROTOCOL)

# Load the graph

eColiNetwork = pickle.load(open('EcoliNetwork.gn', 'r'))


#### PRINT INFOS ON THE NETWORK ########################################


# undirected version of E. coli's network

undirected_net = eColiNetwork.to_undirected()

# Stats

print '-- General informations : \n', \
        nx.info(eColiNetwork)

connected_components = nx.connected_components(undirected_net)

print 'Number of connected components : ', \
        len(connected_components )

print 'Length of the connected components : ', \
        [ len(c) for c in connected_components]

main_component = nx.connected_component_subgraphs(undirected_net)[0]

print "-- Average shortest path in E. coli's main connected component : ", \
       nx.average_shortest_path_length(main_component)

print "-- Diameter of E. coli's main connected component : ", \
       nx.diameter(main_component)

# Infos on a specific gene :

print "-- Genes regulating acs : ", \
       " ".join( eColiNetwork.predecessors('acs') )

print "-- Genes regulated by allr : ", \
       " ".join( eColiNetwork.successors('allr') )


###########  D R A W  ##################################################


from pylab import *

figure()
component = nx.connected_component_subgraphs(undirected_net)[1]
nx.draw(component, node_size = 1000)
savefig('test.png')
#show()
