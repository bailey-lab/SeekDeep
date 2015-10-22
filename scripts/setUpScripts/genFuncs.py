#!/usr/bin/env python

import shutil, os, argparse, sys, stat
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyUtils"))
from utils import Utils

class genHelper:
    @staticmethod
    def generateCompfileFull(outFileName, externalDirLoc, cc, cxx, outName, installDirName, installDirLoc, neededLibs):
        availableLibs = ["CPPITERTOOLS","CPPPROGUTILS","ZI_LIB","BOOST","R","BAMTOOLS","CPPCMS","MATHGL","ARMADILLO",
                         "MLPACK","LIBLINEAR","PEAR","CURL","GTKMM", "BIBSEQ", "BIBCPP", "SEEKDEEP", 
                         "BIBSEQDEV", "SEEKDEEPDEV", "CATCH", "JSONCPP",
                          "TWOBIT", "SEQSERVER","NJHRINSIDE", "PSTREAMS", "MONGOC", "MONGOCXX"]
        neededLibraries = {}
        for lib in neededLibs:
            if ":" in lib:
                libSplit = lib.split(":")
                neededLibraries[libSplit[0].upper()] = libSplit[1]
            else:
                neededLibraries[lib.upper()] = ""
        """
            @todo: Make some of these default to an envirnment CC and CXX and maybe even CXXFLAGS as well 
            @todo: Make availableLibs a more universal constant
        """
        with open(outFileName, "w") as f:
            f.write("CC = {CC}\n".format(CC = cc))
            f.write("CXX = {CXX}\n".format(CXX = cxx))
            f.write("CXXOUTNAME = {NAME_OF_PROGRAM}\n".format(NAME_OF_PROGRAM = outName))
            #f.write("CXXFLAGS = -std=c++11\n")
            f.write("CXXFLAGS = -std=c++14\n")
            f.write("CXXFLAGS += -Wall -ftemplate-depth=1024\n")
            f.write("CXXOPT += -O2 -funroll-loops -DNDEBUG  \n")
            f.write("ifneq ($(shell uname -s),Darwin)\n")
            f.write("\tCXXOPT += -march=native -mtune=native\n" )
            f.write("endif\n")
            f.write("\n")
            f.write("#debug\n")
            f.write("CXXDEBUG = -g -gstabs+ \n")
            f.write("INSTALL_DIR={INSTALL_LOCATION}\n".format(INSTALL_LOCATION = os.path.join(installDirLoc,installDirName)))
            f.write("EXT_PATH=$(realpath {EXTERNAL})\n".format(EXTERNAL = externalDirLoc))
            #f.write("SCRIPTS_DIR=$(realpath scripts)\n")
            f.write("\n")
            for lib in availableLibs:
                if lib in neededLibraries:
                    if neededLibraries[lib] == "":
                        f.write("USE_{LIB} = 1\n".format(LIB = lib))
                    else:
                        f.write("USE_{LIB} = 1#{BRANCH}\n".format(LIB = lib, BRANCH = neededLibraries[lib]))
                else:
                    f.write("USE_{LIB} = 0\n".format(LIB = lib))
                    


    @staticmethod            
    def determineCC(args):
        defaultCC = "clang-3.5"
        if Utils.isMac():
            defaultCC = "clang"
        if not args.CC:
            eCC = os.getenv("CC")
            if(eCC):
                defaultCC =  eCC
        else:
            defaultCC = args.CC[0]
        return defaultCC
    @staticmethod
    def determineCXX(args):
        defaultCXX = "clang++-3.5"
        if Utils.isMac():
            defaultCXX = "clang++"
        if not args.CXX:
            eCXX = os.getenv("CXX")
            if  eCXX:
                defaultCXX = eCXX
        else:
            defaultCXX = args.CXX[0]
        return defaultCXX
    
    @staticmethod
    def parseNjhConfigureArgs():
        parser = argparse.ArgumentParser()
        parser.add_argument('-prefix', type=str, nargs=1)
        parser.add_argument('-externalLibDir', type=str, nargs=1)
        parser.add_argument('-CC', type=str, nargs=1)
        parser.add_argument('-CXX', type=str, nargs=1)
        return parser.parse_args()
    
    @staticmethod
    def mkConfigCmd(name,libs, argv):
        if libs == "":
            cmd = "./scripts/setUpScripts/njhConfigure.py -name {name} ".format(name=name)
        else:
            cmd = "./scripts/setUpScripts/njhConfigure.py -name {name} -libs {libs}".format(name=name, libs=libs)
        if len(sys.argv) > 1:
            cmd += " " + " ".join(sys.argv[1:])
        return cmd


    @staticmethod
    def mkConfigFileStr(name, libs):
        ret = """#!/usr/bin/env python
import shutil, os, argparse, sys, stat
sys.path.append("scripts/pyUtils")
sys.path.append("scripts/setUpScripts")
from utils import Utils
from genFuncs import genHelper

def main():
    name = "{name}"
    libs = "{libs}"
    args = genHelper.parseNjhConfigureArgs()
    cmd = genHelper.mkConfigCmd(name, libs, sys.argv)
    Utils.run(cmd)
    
main()
""".format(name = name, libs = libs)
        return ret
    
    