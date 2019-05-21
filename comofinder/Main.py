#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys
import datetime

from .ReadData import Reader
from .Utils import getArgsList, Fprint
from .Calculate import BPSignificance
from .NetworkRandom import BatchProcess
from .MotifDiscovery import MotifMapReduce


def main_run(argument):
    matx, nodelist, nodemap = Reader(argument.inputfile)
    motif = MotifMapReduce(matx, nodelist)
    motif.motifsize = argument.motifysize
    motif.np = argument.nrofprocess
    motif.run()
    # network randomizing
    Fprint("========Started to generate and process the randomized network========")
    BatchProcess.setSubgfiledir(argument.batchfiledir)
    for i in range(argument.ensemblesize):
        bp = BatchProcess(argument.motifysize, i, argument.nrofprocess)
        bp.randomnetwork(matx, nodelist, nodemap)
    Fprint("========Finished processing the randomzied network========")
    # 统计显著性
    bpcs = BPSignificance(motif.subgraph)
    bpcs.analyzbatch(argument.batchfiledir)
    bpcs.statisticalAnalysis()
    bpcs.outputMotif(argument.motifilename)


def public_run(arg):
    argss = getArgsList(arg)
    main_run(argss)


if __name__ == '__main__':
    st = datetime.datetime.now()
    args = getArgsList(sys.argv[1:])
    main_run(args)
    end = datetime.datetime.now()
    Fprint("The whole algorithm took " + str((end - st).microseconds) + " ms")
