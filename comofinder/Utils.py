# coding=utf-8
import getopt
import os
import sys
import warnings

from .Conf import VERSION, LN, DEBUG
# from Logger import Logs
from .DataStruct import Arguments


def createlist(size, elem=0.0):
    newlsit = []
    for i in range(size):
        newlsit.append(elem)
    return newlsit


def getArgsList(argslist):
    try:
        opts, arg = getopt.getopt(argslist, "hvm:n:f:o:b:e:i:r:", ["version", "help"])
    except getopt.GetoptError as e:
        raise e
    args = Arguments()

    for opt, value in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-v", "--version"):
            print("comofinder python version %s" % VERSION)
            sys.exit(0)
        elif opt == "-m":
            args.motifysize = int(value)
        elif opt == "-n":
            args.nrofprocess = int(value)
        elif opt == "-f":
            args.inputfile = value
        elif opt == "-o":
            args.outputdir = value
        elif opt == "-b":
            args.batchfiledir = value
        elif opt == "-e":
            args.ensemblesize = int(value)
        elif opt == "-i":
            args.iterationsmax = int(value)
        elif opt == "-r":
            args.motifilename = value
    args.checkargs()
    return args


# the method for [stdio(standard input/ouput), log(logging to file), warn(warnning info), error(error)]
# the level for [d(debug), i(info), w(warning), e(error), c(critical)]
def Fprint(msg, method="stdio", level="i"):
    # if method == "stdio":
    #     if level == "i" and DEBUG:
    #         print(msg)
    #     elif level == "d" and DEBUG:
    #         print(msg)
    # elif method == "warn":
    #     warnings.warn(msg)
    # change logger to MotifSystem.web.logger.mslogger
    pass


def IsDir(path, mod="dir"):
    res = False
    if mod == "dir":
        res = os.path.isdir(path)
    elif mod == "file":
        res = os.path.isfile(path)

    return res


def usage():
    if LN == "zh":
        usage_zh()
    elif LN == "en":
        usage_en()


def usage_zh():
    print("\nComoFinder Python %s\n"
          "\t-h, --help:		命令行帮助\n"
          "\t-v, --version:		查看版本\n"
          "\t-o:	指定检测到的图案的输出目录，例如:output/subgraphList\n"
          "\t-m:	指定算法的主题大小，默认值为3\n"
          "\t-n:	指定算法中使用的线程数，默认情况下，程序会检测正在运行的计算机上的可用处理器，并为每个处理器分配子任务\n"
          "\t-f:	指定算法的输入文件，例如“sample_input_network.txt”。必须提供此参数\n"
          "\t-b:	指定随机网络集合的复合子图计数的输出目录，例如:output。必须提供此参数\n"
          "\t-e:	指定随机网络集合的大小，默认值为1000\n"
          "\t-i:	指定每个随机化过程中边缘交换过程中的最大迭代次数，默认值为100\n"
          "\t--mfn:	指定最终的motif文件名，例如“output/motif.txt”。必须提供此参数\n"
          % VERSION
          )


def usage_en():
    print("\nComofinder Python %s\n"
          "\t-m: specify the motif size for the algorithm, e.g. 4. The default value is 3.\n"
          "\t-n: specify the number of threads used in the algorithm,By default,\n"
          "\t\tthe program detects the available processors on the running computer and assigns subtasks to each of them.\n"
          "\t-f: specify the input file for the algorithm, e.g. “sample_input_network.txt”. This parameter must be provided.\n"
          "\t-o: specify the output directory for the detected motifs, e.g. output/subgraphList.\n"
          "\t-b: specify the output directory for the composite subgraphs count of the random networks ensemble,\n"
          "\t\te.g. “output”. This parameter must be provided.\n"
          "\t-e: specify the size of the random network ensemble, e.g. 500. Default value is 1000.\n"
          "\t-i: specify the maximal number of iterations in the edge swap process within each randomization procedure,Default value is 100.\n"
          "\t-mfn: specify the final motif file name, e.g. “output/motif.txt”. This parameter must be provided.\n"
          % VERSION
          )
