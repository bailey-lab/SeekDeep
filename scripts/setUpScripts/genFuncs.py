#!/usr/bin/env python

import shutil, os, argparse, sys, stat
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyUtils"))
from utils import Utils

class genHelper:
    @staticmethod
    def generateCompfileFull(outFileName, externalDirLoc, cc, cxx, outName, installDirName, installDirLoc, neededLibs,ldFlags = "", cxxFlags = ""):
        availableLibs = ["CPPITERTOOLS","CPPPROGUTILS","ZI_LIB","BOOST",
                         "R","BAMTOOLS","CPPCMS","MATHGL","ARMADILLO",
                         "MLPACK","LIBLINEAR","CURL","GTKMM", "BIBSEQ",
                          "BIBCPP", "SEEKDEEP", 
                         "BIBSEQDEV", "SEEKDEEPDEV", "CATCH", "JSONCPP",
                          "TWOBIT", "SEQSERVER","NJHRINSIDE", "PSTREAMS",
                           "MONGOC", "MONGOCXX", "SHAREDMUTEX",
                           "MAGIC", "HTS", "RESTBED", "LIBPCA", "BOOST_FILESYSTEM", 
                           "ELUCIDATOR"]
        neededLibraries = {}
        for lib in neededLibs:
            if ":" in lib:
                libSplit = lib.split(":")
                neededLibraries[libSplit[0].upper()] = libSplit[1]
            else:
                neededLibraries[lib.upper()] = ""

        with open(outFileName, "w") as f:
            f.write("CC = {CC}\n".format(CC = cc))
            f.write("CXX = {CXX}\n".format(CXX = cxx))
            f.write("CXXOUTNAME = {NAME_OF_PROGRAM}\n".format(NAME_OF_PROGRAM = outName))
            f.write("CXXFLAGS = -std=c++14\n")
            f.write("CXXFLAGS += -Wall -ftemplate-depth=1024 -Werror=uninitialized -Werror=return-type -Wno-missing-braces\n")
            if "" != cxxFlags:
                if cxxFlags.startswith("\\"):
                    cxxFlags = cxxFlags[1:]
                f.write("CXXFLAGS += " + cxxFlags +"\n")
            if "" != ldFlags:
                f.write("LD_FLAGS = ")
                if not ldFlags.startswith("-"):
                    f.write("-")
                f.write("{ld_flags}\n".format(ld_flags = " ".join(ldFlags.split(","))))
            f.write("CXXOPT += -O2 -funroll-loops -DNDEBUG  \n")
            f.write("\n")
            f.write("#debug\n")
            f.write("CXXDEBUG = -g -gstabs+ \n")
            f.write("INSTALL_DIR={INSTALL_LOCATION}\n".format(INSTALL_LOCATION = os.path.join(installDirLoc,installDirName)))
            f.write("EXT_PATH=$(realpath {EXTERNAL})\n".format(EXTERNAL = externalDirLoc))
            f.write("\n")
            for lib in availableLibs:
                if lib in neededLibraries:
                    if neededLibraries[lib] == "":
                        f.write("USE_{LIB} = 1\n".format(LIB = lib))
                    else:
                        f.write("USE_{LIB} = 1#{BRANCH}\n".format(LIB = lib, BRANCH = neededLibraries[lib]))
                #else:
                #    f.write("USE_{LIB} = 0\n".format(LIB = lib))
                    


    @staticmethod            
    def determineCC(args, defaultCC = "clang-3.8"):
        if Utils.isMac():
            defaultCC = "clang"
        if not args.CC:
            eCC = os.getenv("CC")
            if(eCC):
                defaultCC =  eCC
        else:
            defaultCC =  Utils.getStrFromStrOrList(args.CC)
        return defaultCC
    
    @staticmethod
    def determineCXX(args, defaultCXX = "clang++-3.8"):
        if Utils.isMac():
            defaultCXX = "clang++"
        if not args.CXX:
            eCXX = os.getenv("CXX")
            if  eCXX:
                defaultCXX = eCXX
        else:
            defaultCXX = Utils.getStrFromStrOrList(args.CXX)
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
    def mkConfigCmd(name,libs, argv, ldflags="", cxxFlags=""):
        if libs == "":
            cmd = "./scripts/setUpScripts/njhConfigure.py -name {name} ".format(name=name)
        else:
            cmd = "./scripts/setUpScripts/njhConfigure.py -name {name} -libs {libs}".format(name=name, libs=libs)
        if "" != ldflags:
            if ldflags.startswith("-"):
                ldflags = ldflags[1:]
            cmd += " -ldFlags " + ldflags
        if "" != cxxFlags:
            addingFlags = " -cxxFlags \""
            if cxxFlags.startswith("-"):
                addingFlags += "\\"
            cmd += addingFlags + cxxFlags + "\""
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
    
    