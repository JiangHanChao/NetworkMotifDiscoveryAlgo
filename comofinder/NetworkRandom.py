# coding=utf-8
import os
import copy
import random

from .Utils import Fprint
from .DataStruct import EdgeNode
from .Conf import OUTPUT_DIR, MAXTRIALS
from .MotifDiscovery import MotifMapReduce


class BatchProcess(object):
    batchdir = None

    def __init__(self, mfsize, index, np):
        self._bfname = str(mfsize) + "_batch_" + str(index) + ".cnt"
        self._mfsize = mfsize
        self._nfp = np

    @staticmethod
    def setSubgfiledir(filedir):
        bfdir = os.path.join(OUTPUT_DIR, filedir)
        if not os.path.exists(bfdir):
            os.makedirs(bfdir)
        elif os.path.isfile(bfdir):
            raise ValueError(bfdir + " not a directory")
        BatchProcess.batchdir = bfdir

    def randomnetwork(self, mat, nlist, nmap):
        currbfpath = os.path.join(BatchProcess.batchdir, self._bfname)
        mdr = MDRmultiedgetype(mat, nmap)
        mdr.run()
        tmdp = MotifMapReduce(mdr.randommatrix, nlist)
        tmdp.motifsize = self._mfsize
        tmdp.np = self._nfp
        # renew the nodeInteraction for each node
        tmdp.renewNodeinteraction()
        tmdp.run()
        tmdp.ouputSubgraph(currbfpath)


class MDRmultiedgetype(object):

    def __init__(self, mat, nodemap):
        self._randomMat = copy.deepcopy(mat)
        self._originMat = mat
        self._nodemap = nodemap

    @property
    def randommatrix(self):
        return self._randomMat

    def run(self):
        eglistmap = self.getdiffedgegroup(self._originMat)
        # shuffle the edges in each group separately
        for key, eglist in eglistmap.items():
            self.shuffEdges(eglist, key.split('_')[1])
        # Further check the network
        self.validatediffegroup()

    def getdiffedgegroup(self, matrix):
        edeglistmap = {}
        rows, cols = matrix.shape
        for i in range(rows):
            for j in range(cols):
                x = matrix[i][j]
                y = matrix[j][i]
                regTypeKey = str(self._nodemap[i]) + str(self._nodemap[j])
                regTypeKey += "_" + str(x) + str(y)
                eg = EdgeNode(i, j)
                # for normal regulations between miRNAs -> genes/TFs and TF -> miRNAs/TFs/genes
                # for x = 1 and y = 0: regTypeKey should be "01_10" or "10_10"
                # for x = 3 and y = 2: regTypeKey should be "11_32" or "12_32"

                # for bidirected regulations between miRNAs/TFs and TFs
                # if x == 1 and y == 1:
                #     if i >= j: continue
                # Skip the self-regulations and only count once for bidirected(undirected) edges
                # regTypeKey should be "01_11" or "11_11"

                # for the gene-gene interactions
                # if x == 2 and y == 2:
                #     if i >= j: continue
                # Skip the self-regulations and only count once for gene interactions
                # regTypeKey should be "22_22"
                if x == y and i >= j:
                    continue

                if regTypeKey in edeglistmap:
                    reglist = edeglistmap[regTypeKey]
                    reglist.append(eg)
                else:
                    reglist = [eg]
                    edeglistmap[regTypeKey] = reglist

        return edeglistmap

    def shuffEdges(self, eglist, egtype):
        failsTrial = 0
        maxdiff = len(eglist)
        shuffeg = []
        while failsTrial < MAXTRIALS and len(eglist) > 0:
            fpos = random.randint(0, len(eglist) - 1)
            swapflag = False
            fedge = eglist[fpos]
            prob = float(len(eglist)) / (len(eglist) + len(shuffeg))
            if random.random() < prob:
                spos = random.randint(0, len(eglist) - 1)
                sedge = eglist[spos]
            else:
                spos = random.randint(0, len(shuffeg) - 1)
                sedge = shuffeg[spos]
                swapflag = True
            # Add one more restriction
            if fedge.src == sedge.src or fedge.des == sedge.des:
                failsTrial += 1
                continue
            cpfedge = EdgeNode(fedge.src, fedge.des)
            cpsedge = EdgeNode(sedge.src, sedge.des)
            # swap the two ends of the two edges
            fedge.des = cpsedge.des
            sedge.des = cpfedge.des
            succflag = True

            if self._randomMat[fedge.src][fedge.des] > 0 or self._randomMat[sedge.src][sedge.des] > 0 \
                    or self._randomMat[fedge.des][fedge.src] > 0 or self._randomMat[sedge.des][sedge.src] > 0 \
                    or fedge.src == fedge.des or sedge.src == sedge.des:
                succflag = False
            elif egtype == "10":
                if self._originMat[fedge.src][fedge.des] > 0 or self._originMat[sedge.src][sedge.des] > 0:
                    succflag = False
            else:
                if self._originMat[fedge.src][fedge.des] > 0 or self._originMat[sedge.src][sedge.des] > 0 \
                        or self._originMat[fedge.des][fedge.src] > 0 or self._originMat[sedge.des][sedge.src] > 0:
                    succflag = False
            if succflag:
                self.updateRandomMatrix(fedge, sedge, cpfedge, cpsedge, egtype)
                if not swapflag:
                    if fpos > spos:
                        eglist.pop(fpos)
                        eglist.pop(spos)
                    else:
                        eglist.pop(spos)
                        eglist.pop(fpos)
                    shuffeg.append(fedge)
                    shuffeg.append(sedge)
                else:
                    eglist.pop(fpos)
                    shuffeg.append(fedge)
                failsTrial = 0
            else:  # Switch the edge back
                fedge.des = cpfedge.des
                sedge.des = cpsedge.des
                failsTrial += 1
        Fprint("Shuffled " + str(maxdiff - len(eglist)) + " number of edges in the current edge list.", level="d")

    def updateRandomMatrix(self, first, second, cpfirst, cpsecond, etype):
        self._randomMat[first.src][first.des] = self._randomMat[cpfirst.src][cpfirst.des]
        self._randomMat[second.src][second.des] = self._randomMat[cpsecond.src][cpsecond.des]
        self._randomMat[cpfirst.src][cpfirst.des] = 0
        self._randomMat[cpsecond.src][cpsecond.des] = 0

        if etype == "10":
            self._randomMat[first.des][first.src] = self._randomMat[cpfirst.des][cpfirst.src]
            self._randomMat[second.des][second.src] = self._randomMat[cpsecond.des][cpsecond.src]
            self._randomMat[cpfirst.des][cpfirst.src] = 0
            self._randomMat[cpsecond.des][cpsecond.src] = 0

    def validatediffegroup(self):
        tempemap = self.getdiffedgegroup(self._randomMat)
        elistmap = self.getdiffedgegroup(self._originMat)
        for key, elist in tempemap.items():
            if len(elist) != len(elistmap[key]):
                raise ValueError("The number of edges in different groups are not the same")
