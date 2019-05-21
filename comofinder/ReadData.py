# coding=utf-8
import os
import numpy as np

from .Conf import FILE_DIR
from .Algorithm import Binarysort
from .DataStruct import EdgeNode, Node


def Reader(filepath):
	fp = os.path.join(FILE_DIR, filepath)
	nodeMap = {}
	edgelist = []
	if not os.path.exists(fp):
		raise ValueError(filepath+" was no found.")
	
	with open(fp, 'r') as f:
		for line in f.readlines():
			infos = line.strip().split('\t')
			oneID = int(infos[0])
			twoID = int(infos[1])
			oneType = int(infos[2])
			twoType = int(infos[3])
			if oneID not in nodeMap:
				nodeMap[oneID] = oneType
			if twoID not in nodeMap:
				nodeMap[twoID] = twoType
			if len(infos) == 4:
				edgeType = 1
			else:
				edgeType = int(infos[4])
			e = EdgeNode(oneID, twoID)
			e.etype = edgeType
			edgelist.append(e)
	return getNodesandEdge(nodeMap, edgelist)


def getNodesandEdge(nodemap, edgelist):
	nodelist = []
	matrix = np.zeros((len(nodemap), len(nodemap)), dtype=np.int)

	for i in range(len(nodemap)):
		atom = Node()
		atom.nodeid = i
		if i not in nodemap:
			raise AttributeError('node id '+str(i)+'no found in node map keys.')
		atom.nodetype = nodemap[i]
		nodelist.append(atom)

	for eg in edgelist:
		matrix[eg.src][eg.des] = eg.etype
		Binarysort(nodelist[eg.src].nodeinteraction, eg.des)
		Binarysort(nodelist[eg.des].nodeinteraction, eg.src)

	return matrix, nodelist, nodemap
