#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import math
import numpy as np

from .Conf import RESULT_DIR, OUTPUT_DIR
from .Utils import IsDir, Fprint, createlist


class BPSignificance(object):

    def __init__(self, subg):
        self._countSubg = subg
        self._pValue = createlist(len(subg.keys()))
        self._zScore = createlist(len(subg.keys()))
        self._ems = None
        self._freqcy = None
        self._meanfreq = None
        self._occupancy = 0
        self._labelindex = {}
        self._labelidxreverse = {}
        self._standeviation = None
        p = 0
        for key in self._countSubg.keys():
            self._labelindex[key] = p
            self._labelidxreverse[p] = key
            p += 1

    def analyzbatch(self, batchdir):
        batchdir = os.path.join(OUTPUT_DIR, batchdir)
        if not IsDir(batchdir):
            raise ValueError(batchdir+" not a directory")
        filelist = os.listdir(batchdir)
        self._ems = len(filelist)
        self._freqcy = np.zeros((self._ems, len(self._countSubg.keys())), dtype=int)

        for i in range(self._ems):
            with open(os.path.join(batchdir, filelist[i]), 'r') as bf:
                bcsubg = {}
                lines = bf.readlines()
                for line in lines:
                    linfo = line.split('\t')
                    bcsubg[linfo[0]] = int(linfo[1])
                for key in self._countSubg.keys():
                    if key in bcsubg:
                        self._freqcy[i][self._labelindex[key]] = bcsubg[key]
                    else:
                        Fprint("Not find " + key + " in the " + filelist[i] + " randomized network", "warn")
        Fprint("Finished read all the batch files")

    def statisticalAnalysis(self):
        for v in self._countSubg.values():
            self._occupancy += v
        self._meanfreq = createlist(len(self._zScore))
        row, col = self._freqcy.shape
        for i in range(row):
            for j in range(col):
                self._meanfreq[j] = self._freqcy[i][j]
                if self._freqcy[i][j] >= self._countSubg[self._labelidxreverse[j]]:
                    self._pValue[j] += 1
        standardDeviation = createlist(len(self._zScore))
        for j in range(col):
            # compute mean frequency
            self._meanfreq[j] /= float(self._ems)
            # compute the ZScore
            Np = float(self._countSubg[self._labelidxreverse[j]])
            for i in range(row):
                mius = self._freqcy[i][j]-self._meanfreq[j]
                standardDeviation[j] += mius*mius
            standardDeviation[j] = math.sqrt(standardDeviation[j]/(row - 1))
            self._zScore[j] = (Np - self._meanfreq[j]) / standardDeviation[j]
            # compute pValue in a more reliable way
            self._pValue[j] = (self._pValue[j] + 1)/(self._ems + 1)
        self._standeviation = standardDeviation
        Fprint("Finish computing the significance")

    def outputMotif(self, filename):
        if not os.path.exists(RESULT_DIR):
            os.makedirs(RESULT_DIR)
        with open(os.path.join(RESULT_DIR, filename), 'w') as f:
            f.write("Subgraph_Label\toccupancy\tmean_frequency\t"
                    "standard_deviation\tZScore\tpValue\t\tfrequency_in_real_network")
            f.write('\n')

            for key, value in self._countSubg.items():
                point = self._labelindex[key]
                f.write(key + "\t")
                f.write(str(value/self._occupancy) + "\t")
                f.write(str(self._meanfreq[point]) + "\t")
                f.write(str(self._standeviation[point]) + "\t")
                f.write(str(self._zScore[point]) + "\t")
                f.write(str(self._pValue[point]) + "\t")
                f.write(str(value))
                f.write("\n")
                point += 1
        Fprint("output motif result wirte finished")
