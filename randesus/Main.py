#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import getopt
import subprocess

from MotifSystem.settings import UPLOAD_RANDESU_ROOT, BASE_DIR


def main_run(argument):
    fpin = os.path.join(UPLOAD_RANDESU_ROOT, 'data', argument.filein)
    fpout = os.path.join(UPLOAD_RANDESU_ROOT, 'result')
    if not os.path.exists(fpout):
        os.makedirs(fpout)
    fpres = os.path.join(fpout, argument.output)
    fnow = os.path.join(BASE_DIR, 'algori', 'randesus', 'randesu')
    cmd = fnow+" "+fpin+" -s "+argument.ms+" -t "+argument.t+" -r "+argument.r
    if argument.d:
        cmd += " -d "
    cmd += " -o "+fpres
    p = subprocess.call(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p != 0:
        raise ValueError("randesu job return code not 0")
    return p


def public_run(args):
    argslist = getArgsList(args)
    main_run(argslist)


class Arguments(object):
    def __init__(self):
        self._filein = None
        self._output = None
        self._ms = None
        self._t = None
        self._r = None
        self._d = None

    @property
    def filein(self):
        return self._filein

    @filein.setter
    def filein(self, value):
        self._filein = value

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, value):
        self._output = value

    @property
    def ms(self):
        return self._ms

    @ms.setter
    def ms(self, value):
        self._ms = value

    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, value):
        self._t = value

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, value):
        self._r = value

    @property
    def d(self):
        return self._d

    @d.setter
    def d(self, value):
        self._d = value


def getArgsList(argslist):
    try:
        opts, arg = getopt.getopt(argslist, "f:r:s:t:o:d:")
    except getopt.GetoptError as e:
        raise e
    args = Arguments()

    for opt, value in opts:
        if opt == "-f":
            args.filein = value
        elif opt == "-r":
            args.r = str(value)
        elif opt == "-s":
            args.ms = str(value)
        elif opt == "-t":
            args.t = str(value)
        elif opt == "-o":
            args.output = value
        elif opt == "-d":
            args.d = value

    return args
