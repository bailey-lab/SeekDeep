#!/usr/bin/env python

import shutil, os, argparse, sys, stat, time
from genFuncs import genHelper
from utils import Utils
from color_text import ColorText as CT
    
def genCommonIncludes(filename, includes):
     with open(filename, "w") as f:
         f.write("#pragma once\n")
         f.write("\n")
         f.write("// Created on "+ time.strftime("%Y/%m/%d") +"\n")
         f.write("// Add system libraries here\n")
         f.write("\n")
         for i in includes:
             f.write("#include <" + i + ">\n")
            
def genTypeDefs(filename, projName):
     with open(filename, "w") as f:
         f.write("#pragma once\n")
         f.write("\n")
         f.write("// Created on "+ time.strftime("%Y/%m/%d") +"\n")
         f.write("// Add project typedefs here\n")
         f.write("\n")
         for i in ["vector", "map", "string"]:
             f.write("#include <" + i + ">\n")
         f.write("\n")
         f.write("namespace " + projName+" {\n")
         f.write("//typedef std::map<std::string, std::string> MapStrStr;\n")
         f.write("//typedef std::vector<std::string> VecStr;\n")
         f.write("}  // namespace " + projName + "\n")
        
def genCommon(dest, projName, includes):
    genCommonIncludes(os.path.join(dest, "src/" + projName + "/common/allSystemIncludes.h"), includes)
    genTypeDefs(os.path.join(dest, "src/" + projName + "/common/typedefs.h"), projName)
    with open(os.path.join(dest, "src/" + projName + "/common.h"), "w") as f:
        f.write("#pragma once\n")
        f.write("\n")
        f.write("// Created on "+ time.strftime("%Y/%m/%d") +"\n")
        f.write("// Including headers in common\n")
        f.write("\n")
        f.write("#include \""+ projName + "/common/allSystemIncludes.h\"\n")
        f.write("#include \""+ projName + "/common/typedefs.h\"\n")
        
def genWholeProjInclude(dest, projName):
     with open(os.path.join(dest, "src/" + projName + ".h"), "w") as f:
         f.write("#pragma once\n")
         f.write("\n")
         f.write("// Created on "+ time.strftime("%Y/%m/%d") +"\n")
         f.write("// Including whole project\n")
         f.write("\n")
         f.write("#include \""+ projName + "/common.h\"")

def genMain(dest, projName):
    with open(os.path.join(dest, "src/" + "main.cpp"), "w") as f:
        f.write("\n")
        f.write("// Created on "+ time.strftime("%Y/%m/%d") +"\n")
        f.write("// main.cpp\n")
        f.write("\n")
        f.write("#include \""+ projName + ".h\"")
        f.write("\n")
        f.write("int main(int argc, char* argv[]){\n")
        f.write("    std::cout << \"Hello " +  projName +"!\" << std::endl;\n")
        f.write("    return 0;\n")
        f.write("}\n")
        
def genSrc(dest, projName, includes):
    os.mkdir(os.path.join(dest, "src"))
    os.mkdir(os.path.join(dest, "src/" + projName))
    os.mkdir(os.path.join(dest, "src/" + projName + "/common"))     
    genCommon(dest, projName,  includes)
    genWholeProjInclude(dest, projName)
    genMain(dest, projName)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-projName', type=str, nargs=1, required=True)
    parser.add_argument('-dest', type=str, nargs=1, required=True)
    parser.add_argument('-CC', type=str, nargs=1)
    parser.add_argument('-CXX', type=str, nargs=1)
    parser.add_argument('-externalLoc', type=str, nargs=1)
    parser.add_argument('-neededLibs', type=str, nargs=1)
    return parser.parse_args()

def main():
    args = parse_args()
    projectOut = os.path.join(args.dest[0], args.projName[0])
    os.mkdir(projectOut)
    genSrc(projectOut, args.projName[0], ["iostream", "string", "unistd.h", "vector", "stdint.h", "stdio.h", "cstddef", "utility", "map", "algorithm"])
    CC = genHelper.determineCC(args)
    CXX = genHelper.determineCXX(args)
    external = "external"
    outname = args.projName[0]
    prefix = "./"
    installName = args.projName[0]
    neededLibs = "none"        
    if args.externalLoc:
        external = os.path.realpath(args.externalLoc[0])
    if args.neededLibs:
        neededLibs = args.neededLibs[0].split(",")
    genHelper.generateCompfileFull(os.path.join(projectOut, "compfile.mk"), external, CC, CXX, outname, installName, prefix, neededLibs)
    
    exFrom = os.path.abspath(os.path.dirname(__file__))
    cpSetUpCmd = exFrom + "/copySetUpFiles.py -from " + exFrom +"/../ -to " + projectOut
    print CT.boldBlack(cpSetUpCmd)
    Utils.run(cpSetUpCmd)
    cpMakefilesCmd = "cp " + exFrom + "/cppSetUpFiles/*akefile* " + projectOut
    print CT.boldBlack(cpMakefilesCmd)
    Utils.run(cpMakefilesCmd)
    
main()
