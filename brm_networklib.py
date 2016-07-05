import sys
import os
import re
import collections as cl
import networkx as nx
import networkxgmml as nxgml

def matrix2network(inMatrix,directed=False,cutoff=-1.0):
	"""Converts a two dimensional matrix into a networkx network object. Edges with weight
below the cutoff will be ignored. Currently only works for undirected networks"""

	oFile = open(inMatrix, 'r')
	f = oFile.readlines()

	#Catches the column labels
	header = f.pop(0)
	header = header.rstrip('\n')
	cNames = list(header.split("\t"))[1:]

	NODES = {}
	EDGES = {}

	#Read data from the matrix
	count = 0
	printed = {}
	if directed == True:
		invComp = ''
	for ln in f:
		ln = ln.rstrip('\n')
		row = ln.split()
		rowName = row.pop(0)
		if len(row) > 0:
			if 'NA' not in row[0]:
				for idx, entry in enumerate(row):
					if float(entry) >= cutoff:
						comp = "{0}.{1}".format(rowName,cNames[idx])
						if directed == False:
							invComp = "{0}.{1}".format(cNames[idx],rowName)

							if rowName != cNames[idx] and not printed.has_key(invComp):
								NODES[rowName] = ''
								NODES[cNames[idx]] = ''

								if not EDGES.has_key(rowName):
									EDGES[rowName] = {}
								EDGES[rowName][cNames[idx]] = float(entry)
								printed[comp] = ''

	#Generate networkx object
	network = nx.Graph()
	for i in NODES:
		network.add_node(i,graphics={"fill":"#FFFED9","w":35,"h":35,"d":35,"type":"circle"})
	for i in EDGES:
		for j in EDGES[i]:
			try:
				network.add_edge(i,j,weight=EDGES[i][j])
			except KeyError:
				print 'could not find edge: {0} -> {1}'.format(i,j)
	return network


def read_cytoscape_location(network,cyt):
	"""Pull x/y coordinates from a cytoscape xml output because it doesn't output things correctly"""

	NODE = {}
	with open(cyt,'r') as f:
		for ln in f:
			if ' <node id=' in ln:
				nodeSearch = re.search("label=\"(\S+)\"",ln)
				if nodeSearch:
					node = nodeSearch.group(1)
					NODE[node] = {"x":"","y":""}

				else:
					print "Couldn't find node in {0}".format(ln)
			elif "<graphics" in ln and "x=" in ln and "y=" in ln:
				m = re.search("x=\"(\S+)\"",ln)
				if m:
					NODE[node]["x"] = m.group(1)
				else:
					print "Couldn't read x coordinate data in {0}".format(ln)
				m = re.search("y=\"(\S+)\"",ln)
				if m:
					NODE[node]["y"] = m.group(1)
				else:
					print "Couldn't read y coordinate data in {0}".format(ln)

	for n in NODE:
		if network.node.has_key(n):
			if NODE[n]["x"] == "" or NODE[n]["y"] == "":
				print "Incomplete coordinate datat for {0}".format(n)
			else:
				network.node[n]['x'] = float(NODE[n]['x'])
				network.node[n]['y'] = float(NODE[n]['y'])

	#Delete nodes that don't have location data
	remove = []
	for idTag in network.node:
		if not network.node[idTag].has_key('x'):
			remove.append(idTag)
		elif network.node[idTag]['x'] == 0 and network.node[idTag]['x'] == 0:
			remove.append(idTag)
	for n in remove:
		network.remove_node(n)
	return network

def currate_edges(network,minVal,maxVal):
	"""Remove edges if their weight is outside a given range"""
	try:
		minVal = float(minVal)
		minVal = float(maxVal)
	except ValueError:
		return "Minimum and maximum values must be numbers"
	for e in network.edges(data=True):
		if e[2]['weight'] < minVal or e[2]['weight'] > maxVal:
			network.remove_edge(e[0],e[1])
	return network


def add_node_attribute(network,inFile,attributeName):
	"""Add attributes to nodes from a tab delimited file"""
	with open(inFile,'r') as f:
		for ln in f:
			ln = ln.replace('\n','')
			ln = ln.split()
			if network.has_node(ln[0]):
				network.node[ln[0]][attributeName] = ln[1]
	return network


def map_node_colors(network,attribute,colorFile):
	"""Colors nodes based on their attributes with colors defined by a tab-delimited file"""
	COLOR = {}
	with open(colorFile,'r') as f:
		for ln in f:
			ln = ln.replace('\n','')
			ln = ln.split()
			COLOR[ln[0]] = ln[1]
	for i in network.nodes():
		if network.node[i].has_key(attribute):
			att = network.node[i][attribute]
			if COLOR.has_key(att):
				network.node[i]['graphics']["fill"] = COLOR[att]
		else:
			network.node[i]['graphics']["fill"] = "#FFFED9"
	return network


def fet_attribute_connectedness(network,att,att2='',attValue='',attValue2='',LIMIT={},minWeight=0.0):
	"""
	Use Fisher's exact test to look for enrichment of attributes and edges between nodes

	Variants:
		1.  Test for a correlation between a shared attribute value and an edge between nodes
		(ie nodes linked to other nodes of the same color)
		Inputs: network, attribute

		2.  Test for a correlation between shared attribute value and an edge between nodes with a 
		specific attribute value (ie specifically blue nodes linked to other blue nodes)
		Inputs: network, attribute, attribute value

		3.  Test for a correlation between edges between nodes with two different values for the
		same attribute (ie blue nodes linked to red nodes)
		Inputs: network, attribute, attribute value, second attribute value

		4.  Test for a correlation between edges between nodes and specific values for two 
		different attributes (ie blue nodes linked to square nodes)
		Inputs: network, attribute, second attribute, attribute value, second attribute value

		#a.  Limit tested nodes to those with particular attribute values
		Additional inputs: Dict of attribute:value pairs
	"""

	from scipy import stats

	#A = match attribute criteria
	#B = edge between nodes
	aANDb = 0
	aNOTb = 0
	bNOTa = 0
	notAB = 0

	if not attValue and not att2:
		if LIMIT:
			print "Running variant 1a"
		else:
			print "Running variant 1"
	elif attValue and not attValue2:
		if LIMIT:
			print "Running variant 2a"
		else:
			print "Running variant 2"
	elif attValue and attValue2 and not att2:
		if LIMIT:
			print "Running variant 3a"
		else:
			print "Running variant 3"
	elif attValue and attValue2 and att2:
		if LIMIT:
			print "Running variant 4a"
		else:
			print "Running variant 4"

	#Fetching required nodes
	nodes = []
	for i in network.nodes():
		if LIMIT:
			keep = 1
			for lim in LIMIT:
				if network.node[i].has_key(lim):
					if not network.node[i][lim] == LIMIT[lim]:
						keep = 0
			if keep == 1:
				nodes.append(i)
		else:
			nodes.append(i)

	#Variant 1
	if not attValue and not att2:
		for idx, i in enumerate(nodes):
			for j in nodes[idx+1:]:
				if network.node[i].has_key(att) and network.node[j].has_key(att):
					if network[i].has_key(j) and network.node[i][att] == network.node[j][att]:
						if network[i][j]['weight'] >= minWeight:
							aANDb += 1
					elif network.node[i][att] == network.node[j][att]:
						aNOTb += 1
					elif network[i].has_key(j):
						if network[i][j]['weight'] >= minWeight:
							bNOTa += 1
					else:
						notAB += 1

	#Variant 2
	elif attValue and not attValue2:
		for idx, i in enumerate(nodes):
			for j in nodes[idx+1:]:
				if network.node[i].has_key(att) and network.node[j].has_key(att):
					if network.node[i][att] == attValue:
						if network[i].has_key(j) and network.node[j][att] == attValue:
							if network[i][j]['weight'] >= minWeight:
								aANDb += 1
						elif network.node[j][att] == attValue:
							aNOTb += 1
						elif network[i].has_key(j):
							if network[i][j]['weight'] >= minWeight:
								bNOTa += 1
						else:
							notAB += 1

	#Variant 3
	elif attValue and attValue2 and not att2:
		for idx, i in enumerate(nodes):
			for j in nodes[idx+1:]:
				if network.node[i].has_key(att) and network.node[j].has_key(att):
					if network.node[i][att] == attValue:
						if network[i].has_key(j) and network.node[j][att] == attValue2:
							if network[i][j]['weight'] >= minWeight:
								aANDb += 1
						elif network.node[j][att] == attValue2:
							aNOTb += 1
						elif network[i].has_key(j):
							if network[i][j]['weight'] >= minWeight:
								bNOTa += 1
						else:
							notAB += 1

	#Variant 4
	elif attValue and attValue2 and att2:
		for idx, i in enumerate(nodes):
			for j in nodes[idx+1:]:
				if network.node[i].has_key(att) and network.node[j].has_key(att2):
					if network.node[i][att] == attValue:
						if network[i].has_key(j) and network.node[j][att2] == attValue2:
							if network[i][j]['weight'] >= minWeight:
								aANDb += 1
						elif network.node[j][att2] == attValue2:
							aNOTb += 1
						elif network[i].has_key(j):
							if network[i][j]['weight'] >= minWeight:
								bNOTa += 1
						else:
							notAB += 1

	oddsRatio, pValue = stats.fisher_exact([[aANDb, bNOTa], [aNOTb, notAB]])
	return oddsRatio,pValue


