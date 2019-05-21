#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import copy
import numpy as np

from .Utils import Fprint
from .Conf import MOTIFTYPE
from .Algorithm import Binarysort, Binarycontain, Getlabel, Swap


class CMFinder(object):

    def __init__(self, matrix, thold, msize):
        self._mat = matrix
        self._nlist = None
        self._subgraph = {}
        self._THRESHOLD = thold
        self._MOTIFSIZE = msize

    def Finder(self, mRNAlist, nodelist):
        Fprint("Start to extend the seed graphs from 0 to %d" % len(mRNAlist))
        self._nlist = nodelist
        for mRNAid in mRNAlist:
            mRNA = nodelist[mRNAid]
            Vext = []
            Vsub = []
            Vopn = []
            for part in mRNA.nodeinteraction:
                if part > mRNAid:
                    Vext.append(part)
                    Vopn.append(part)
            Vsub.append(mRNAid)
            self.extendGsub(Vsub, Vext, Vopn, mRNAid)
            del Vext, Vsub, Vopn
        Fprint("Finish extending the current seed graphs...")

        return self._subgraph

    def extendGsub(self, Vsub, Vext, Vopn, V):
        Fprint("vsub:", level="d")
        Fprint(Vsub, level="d")
        Fprint("vext:", level="d")
        Fprint(Vext, level="d")
        Fprint("vopen:", level="d")
        Fprint(Vopn, level="d")

        while len(Vext) != 0:
            W = Vext.pop(0)
            flag = -1
            if len(Vsub) == self._MOTIFSIZE - 2:
                mark = [False, False, False]
                if len(Vsub) > 1:
                    for v in Vsub:
                        ntype = self._nlist[v].nodetype
                        mark[ntype] = True

                ntype = self._nlist[W].nodetype
                mark[ntype] = True
                if not mark[1]:
                    flag = 1
                elif not mark[2]:
                    flag = 2

            # find all the nodes u from Nexcl(w,Vsubgraph) and u > v
            newVext = copy.copy(Vext)
            newVopn = copy.copy(Vopn)
            self.findInNexcl(W, V, Vsub, newVext, newVopn, flag)

            if len(newVext) == 0:
                continue
            newVsub = copy.copy(Vsub)
            Binarysort(newVsub, W)

            if len(newVsub) == self._MOTIFSIZE - 1:
                self.countGsub(newVsub, newVext)
            else:
                self.extendGsub(newVsub, newVext, newVopn, V)

    def findInNexcl(self, w, v, Vsub, Vext, Vopn, flag):
        # remove those nodes that do not meet our requirements from the Vext
        if flag != -1:
            Vetcp = copy.copy(Vext)
            for vt in Vetcp:
                if self._nlist[vt].nodetype != flag:
                    Vext.remove(vt)
            del Vetcp

        for part in self._nlist[w].nodeinteraction:
            if part > v and Binarycontain(Vsub, part) < 0 and not Binarysort(Vopn, part):
                if flag == -1 or self._nlist[part].nodetype == flag:
                    Binarysort(Vext, part)

    def countGsub(self, Vsub, Vext):
        for ve in Vext:
            newVsub = copy.copy(Vsub)
            newVsub.append(ve)
            self.countGSubs(newVsub)

    def countGSubs(self, Gsub):
        # SubGraph Canonical Label
        Gsublabel = Getlabel(self._mat, self._nlist, Gsub)
        if Gsublabel in self._subgraph:
            count = self._subgraph[Gsublabel]
            count += 1
            self._subgraph[Gsublabel] = count
        else:
            self._subgraph[Gsublabel] = 1


class Rmvisomorphism(object):

    def __init__(self, subG, thold):
        self._threshold = thold
        self._subgraph = subG
        self._newsubg = {}

    def Remove(self):
        for label, value in self._subgraph.items():
            canonicalabel = self.getCanonicalabel(label)
            if canonicalabel in self._newsubg:
                newcnt = self._newsubg[canonicalabel]
                value += newcnt
            self._newsubg[canonicalabel] = value

        return self._newsubg

    @staticmethod
    def getCanonicalabel(sublabel):
        matlabel, typelabel = sublabel.split('_')
        tlen = len(typelabel)
        matrix = np.zeros((tlen, tlen), dtype=np.int)
        row, col = matrix.shape
        for i in range(row):
            for j in range(col):
                matrix[i][j] = int(matlabel[i * tlen + j])

        # type list array
        typegroup = []
        typeComb = []
        for i in range(MOTIFTYPE):
            typegroup.append([])
        for i in range(tlen):
            typegroup[int(typelabel[i])].append(i)
        # Get the canonical labels: maximum lexicographic string
        for group in typegroup:
            pm = Permutate()
            comb = pm.getCombination(group)
            typeComb.append(comb)

        em = Enumerate()
        shuffle = em.getEnumeration(typeComb)
        maxlabel = ""
        for shff in shuffle:
            currlabel = ""
            nodetype = ""
            for i in range(len(shff)):
                for j in range(len(shff)):
                    currlabel += str(matrix[shff[i]][shff[j]])
                nodetype += str(typelabel[shff[i]])
            currlabel += "_" + nodetype
            if currlabel > maxlabel:
                maxlabel = currlabel
        return maxlabel


class Permutate(object):

    def __init__(self):
        self._comba = []

    def getCombination(self, subg):
        self.arrange(subg, 0, len(subg))
        return self._comba

    def arrange(self, substr, start, lens):
        if start == lens - 1:
            newstr = copy.copy(substr)
            self._comba.append(newstr)
        else:
            for i in range(start, lens):
                Swap(substr, start, i)
                self.arrange(substr, start + 1, lens)
                Swap(substr, start, i)


class Enumerate(object):
    def __init__(self):
        self._enumlist = []

    def getEnumeration(self, typecomb):
        altogether = []
        self.iterative(0, typecomb, altogether)
        return self._enumlist

    def iterative(self, point, combina, altogether):
        if point == len(combina):
            self._enumlist.append(altogether)
            return
        subset = combina[point]
        for ss in subset:
            cp = self.addValue(ss, altogether)
            self.iterative(point + 1, combina, cp)

    @staticmethod
    def addValue(subset, altogether):
        cp = copy.copy(altogether)
        for s in subset:
            cp.append(s)
        return cp
