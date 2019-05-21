#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import getopt
import subprocess

from MotifSystem.settings import UPLOAD_WARSWAP_ROOT, BASE_DIR


def main_run(argument):
    args0 = os.path.join(UPLOAD_WARSWAP_ROOT, 'data', argument.edgef)
    args1 = os.path.join(UPLOAD_WARSWAP_ROOT, 'data', argument.vtxf)
    args2 = argument.motifsize
    args3 = os.path.join(UPLOAD_WARSWAP_ROOT, 'result')
    if not os.path.exists(args3):
        os.makedirs(args3)
    args3 = os.path.join(args3, argument.output)
    args4 = os.path.join(UPLOAD_WARSWAP_ROOT, 'output')
    args5 = argument.tmp
    fnow = os.path.join(BASE_DIR, 'algori', 'warswap', 'cmdmain.jar')
    cmd = "java -jar " + fnow + " "+args0+" "+args1+" "+args2+" "+args3+" "+args4+" "+args5
    p = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p != 0:
        raise ValueError("warswap job return code not 0")
    return p


def public_run(args):
    argslist = getArgsList(args)
    main_run(argslist)


class Arguments(object):
    def __init__(self):
        self._edgef = None
        self._vtxf = None
        self._motifsize = None
        self._output = None
        self._tmp = None

    @property
    def edgef(self):
        return self._edgef

    @edgef.setter
    def edgef(self, value):
        self._edgef = value

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, value):
        self._output = value

    @property
    def vtxf(self):
        return self._vtxf

    @vtxf.setter
    def vtxf(self, value):
        self._vtxf = value

    @property
    def motifsize(self):
        return self._motifsize

    @motifsize.setter
    def motifsize(self, value):
        self._motifsize = value

    @property
    def tmp(self):
        return self._tmp

    @tmp.setter
    def tmp(self, value):
        self._tmp = value


def getArgsList(argslist):
    try:
        opts, arg = getopt.getopt(argslist, "e:v:m:o:t:")
    except getopt.GetoptError as e:
        raise e
    args = Arguments()

    for opt, value in opts:
        if opt == "-e":
            args.edgef = value
        elif opt == "-v":
            args.vtxf = value
        elif opt == "-m":
            args.motifsize = str(value)
        elif opt == "-o":
            args.output = value
        elif opt == "-t":
            args.tmp = value

    return args
