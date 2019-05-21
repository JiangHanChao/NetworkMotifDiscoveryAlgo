# -*- coding: utf-8 -*-


class Arguments(object):
    """Arguments: cmd run args list"""

    def __init__(self):
        self._motifysize = 3
        self._nrofprocess = -1
        self._ensemblesize = 500
        self._iterationsmax = 100
        self._motifilename = "motify_result.txt"
        self._inputfile = None
        self._outputdir = None
        self._batchfiledir = None

    # ms
    @property
    def motifysize(self):
        return self._motifysize

    @motifysize.setter
    def motifysize(self, value):
        self._motifysize = value

    # nr
    @property
    def nrofprocess(self):
        return self._nrofprocess

    @nrofprocess.setter
    def nrofprocess(self, value):
        self._nrofprocess = value

    # fp
    @property
    def inputfile(self):
        return self._inputfile

    @inputfile.setter
    def inputfile(self, value):
        self._inputfile = value

    # opd
    @property
    def outputdir(self):
        return self._outputdir

    @outputdir.setter
    def outputdir(self, value):
        self._outputdir = value

    # bfd
    @property
    def batchfiledir(self):
        return self._batchfiledir

    @batchfiledir.setter
    def batchfiledir(self, value):
        self._batchfiledir = value

    # ens
    @property
    def ensemblesize(self):
        return self._ensemblesize

    @ensemblesize.setter
    def ensemblesize(self, value):
        self._ensemblesize = value

    # ite
    @property
    def iterationsmax(self):
        return self._iterationsmax

    @iterationsmax.setter
    def iterationsmax(self, value):
        self._iterationsmax = value

    # mfn
    @property
    def motifilename(self):
        return self._motifilename

    @motifilename.setter
    def motifilename(self, value):
        self._motifilename = value

    def __str__(self):
        return 'motif size:%d\tnrofprocess:%d\tensemble size:%d\tmax iterations:%d\nmotif filename:%s\t' \
               'input file:%s\toutput dir:%s\tbatchfile dir:%s' % (
                   self._motifysize, self._nrofprocess, self._ensemblesize, self._iterationsmax,
                   self._motifilename, self._inputfile, self._outputdir, self._batchfiledir)

    # public function
    def checkargs(self):
        if None in (self._inputfile, self._batchfiledir):
            raise NameError("some parameter must be provided, -h for help")


class EdgeNode(object):
    def __init__(self, src, des):
        self._des = des
        self._src = src
        self._etype = -1

    @property
    def des(self):
        return self._des

    @des.setter
    def des(self, value):
        self._des = value

    @property
    def src(self):
        return self._src

    @src.setter
    def src(self, value):
        self._src = value

    @property
    def etype(self):
        return self._etype

    @etype.setter
    def etype(self, value):
        self._etype = value

    def __str__(self):
        return 'des:%d\tsrc:%d\ttype:%d' % (self._des, self._src, self._etype)


class Node(object):
    """composite  node"""

    def __init__(self):
        self._nodename = None
        self._nodeid = None
        self._indg = None
        self._outdg = None
        self._actualid = None
        self._degree = None
        self._nodeinteraction = []
        self._geneinteraction = []
        self._nodetype = -1
        self._edgetype = -1
        self._weight = 0.0
        self._capacity = None
        self._sampleweight = 0.0
        self._familyid = None

    @property
    def nodename(self):
        return self._nodename

    @nodename.setter
    def nodename(self, value):
        self._nodename = value

    # node index
    @property
    def nodeid(self):
        return self._nodeid

    @nodeid.setter
    def nodeid(self, value):
        self._nodeid = value

    # node in-degree
    @property
    def indg(self):
        return self._indg

    @indg.setter
    def indg(self, value):
        self._indg = value

    @property
    def outdg(self):
        return self._outdg

    @outdg.setter
    def outdg(self, value):
        self._outdg = value

    # This variant stores the nodeID before filtering
    @property
    def actualid(self):
        return self._actualid

    @actualid.setter
    def actualid(self, value):
        self._actualid = value

    @property
    def degree(self):
        return self._degree

    @degree.setter
    def degree(self, value):
        self._degree = value

    # this is for regulations
    @property
    def nodeinteraction(self):
        return self._nodeinteraction

    @nodeinteraction.setter
    def nodeinteraction(self, value):
        self._nodeinteraction.append(value)

    # this is specifically for gene interactions
    @property
    def geneinteraction(self):
        return self._geneinteraction

    @geneinteraction.setter
    def geneinteraction(self, value):
        self._geneinteraction.append(value)

    # 0: microRNA; 1: transcription factor; 2: target genes
    @property
    def nodetype(self):
        return self._nodetype

    @nodetype.setter
    def nodetype(self, value):
        self._nodetype = value

    # 0: no edges; 1: regulations; 2: interactions; 3: regulations plus interactions
    @property
    def edgetype(self):
        return self._edgetype

    @edgetype.setter
    def edgetype(self, value):
        self._edgetype = value

    @property
    def weight(self):
        return self._weight

    @weight.setter
    def weight(self, value):
        self._weight = value

    # The following two members are for the WaRSwap randomization algorithm
    @property
    def capacity(self):
        return self._capacity

    @capacity.setter
    def capacity(self, value):
        self._capacity = value

    @property
    def sampleweight(self):
        return self._sampleweight

    @sampleweight.setter
    def sampleweight(self, value):
        self._sampleweight = value

    @property
    def familyid(self):
        return self._familyid

    @familyid.setter
    def familyid(self, value):
        self._familyid = value

    def __str__(self):
        return 'node name:%s\tnode id:%d\tnode type:%d\tnode weight:%d' % (
            self._nodename, self._nodeid, self._nodetype, self._weight)
