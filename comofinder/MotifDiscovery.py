#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import math

from .Utils import Fprint
from datetime import datetime
from .Algorithm import Binarysort
from multiprocessing import cpu_count
from .Comofinder import CMFinder, Rmvisomorphism


class MotifMapReduce(object):
    def __init__(self, matrix, nodelist):
        self._subgraph = None
        self._matrix = matrix
        self._nodelist = nodelist
        self._mRNAlist = []
        self._TFlist = []
        self._genelist = []
        self._motifsize = None
        self._np = None

        for atom in nodelist:
            if atom.nodetype == 0:
                self._mRNAlist.append(atom.nodeid)
            elif atom.nodetype == 1:
                self._TFlist.append(atom.nodeid)
            elif atom.nodetype == 2:
                self._genelist.append(atom.nodeid)
            else:
                Fprint("the " + str(atom.nodeid) + "th node type no found, skip", "warn")

    @property
    def subgraph(self):
        return self._subgraph

    @property
    def motifsize(self):
        return self._motifsize

    @motifsize.setter
    def motifsize(self, value):
        self._motifsize = value

    @property
    def np(self):
        return self._np

    @np.setter
    def np(self, value):
        self._np = value

    def run(self):
        if self.np < 0:
            self._np = cpu_count()
        assert self.np > 0
        Fprint("processors gets " + str(self.np), level="d")

        threshold = math.ceil(len(self._mRNAlist) / self.np)
        Fprint(("THRESHOLD is %d" % threshold), level="d")
        # enumerate subgrahp network
        tstart = datetime.now()
        finder = CMFinder(self._matrix, threshold, self.motifsize)
        GsubReslut = finder.Finder(self._mRNAlist, self._nodelist)
        tend = datetime.now()
        Fprint("The enumeration subgraph took " + str((tend - tstart).microseconds) + " ms")

        # reomve isomorphism subgrahp
        tstart = datetime.now()
        threshold = math.ceil(len(GsubReslut) / self.np)
        if threshold < self.np:
            threshold = len(GsubReslut)
        Fprint(("Remove THRESHOLD is %d" % threshold), level="d")

        rmfinder = Rmvisomorphism(GsubReslut, threshold)
        self._subgraph = rmfinder.Remove()
        del GsubReslut
        tend = datetime.now()
        Fprint("The isomorphism removal took " + str((tend - tstart).microseconds) + " ms")

    def renewNodeinteraction(self):
        for node in self._nodelist:
            node.nodeinteraction.clear()
        for src in self._nodelist:
            for des in self._nodelist:
                if self._matrix[src.nodeid][des.nodeid] != 0:
                    Binarysort(src.nodeinteraction, des.nodeid)
                    Binarysort(des.nodeinteraction, src.nodeid)

    def ouputSubgraph(self, filename):
        with open(filename, "w") as f:
            for key, value in self._subgraph.items():
                f.write(key + "\t" + str(value))
                f.write('\n')
