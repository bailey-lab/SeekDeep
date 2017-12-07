#!/usr/bin/env python2

import shutil, os, argparse, sys, stat, time
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "setUpScripts"))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyUtils"))
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

def genReadme(projDir, projName, overwrite = False):
    readmeFilename = os.path.join(projDir, "README.md");
    if os.path.exists(readmeFilename) and not overwrite:
        raise Exception("File " + readmeFilename + "already exists")
    with open(readmeFilename, "w") as f:
        f.write("#" + projName + "\n")
        f.write("run the following to configure and download libraries\n\n")
        f.write("```bash\n")
        f.write("./configure.py\n")
        f.write("./setup.py --compfile compfile.mk --outMakefile makefile-common.mk\n")
        f.write("```\n")
        f.write("#Use already existing External directory\n\n")
        f.write("```bash\n")
        f.write("./configure.py -externalLibDir [EXTERNAL_DIR]\n")
        f.write("./setup.py --compfile compfile.mk --outMakefile makefile-common.mk\n")
        f.write("```\n")
        

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--projName', type=str, required=True, help = "Name of the project, a directory will be created with this name")
    parser.add_argument('--dest', type=str, required=True, help = "Destionation of the porject direcotry")
    parser.add_argument('--namespace', type=str, help = "Name of the namespace for the project, otherwise assumes the project name")
    parser.add_argument('--overwrite', action = "store_true", help = "Overwrite if project directory already exists")
    parser.add_argument('--CC', type=str, help = "c compiler")
    parser.add_argument('--CXX', type=str, help = "c++ compiler")
    parser.add_argument('--externalLoc', type=str, help = "A location of a external lib where cpp libraries are already downlaoded")
    parser.add_argument('--neededLibs', type=str, help = "Libraries you want to install, e.g. boost:1_60_0")
    
    return parser.parse_args()

def genCppProject(args):
    projectOut = os.path.join(args.dest, args.projName)
    if os.path.exists(projectOut):
        if args.overwrite:
            shutil.rmtree(projectOut)
        else:
            raise Exception("Directory " + str(projectOut) + " already exists, use --overWrite to delete")
    #create project dir
    os.mkdir(projectOut)
    #generate skeleton source code directory
    genSrc(projectOut, args.projName, ["iostream", "string", "unistd.h", "vector", "cstdint", "cstdio", "cstddef", "utility", "map", "unordered_map", "algorithm"])
    #determine c++ and c compilers
    CC = genHelper.determineCC(args)
    CXX = genHelper.determineCXX(args)
    external = "external"
    outname = args.projName
    prefix = "./"
    installName = args.projName
    neededLibs = "none"        
    if args.externalLoc:
        external = os.path.realpath(args.externalLoc)
    if args.neededLibs:
        neededLibs = args.neededLibs.split(",")
    #generate the compfile
    genHelper.generateCompfileFull(os.path.join(projectOut, "compfile.mk"), external, CC, CXX, outname, installName, prefix, neededLibs)
    #generate config file
    with open(os.path.join(projectOut, "configure.py"), "w") as configFile:
        if args.neededLibs:
            configFile.write(genHelper.mkConfigFileStr(outname, args.neededLibs))
        else:
            configFile.write(genHelper.mkConfigFileStr(outname, ""))
    #make executable
    os.chmod(os.path.join(projectOut, "configure.py"), stat.S_IXGRP | stat.S_IXOTH | stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IWUSR)
    #determine this file's location to dtermine where to copy setup and make files to
    exFrom = os.path.abspath(os.path.dirname(__file__))
    cpSetUpCmd = exFrom + "/copySetUpFiles.py -from " + exFrom + "/../../ -to " + projectOut
    print CT.boldBlack(cpSetUpCmd)
    Utils.run(cpSetUpCmd)
    cpMakefilesCmd = "cp " + exFrom + "/../cppMakefiles/Makefile " + projectOut
    print CT.boldBlack(cpMakefilesCmd)
    Utils.run(cpMakefilesCmd)
    #generate README.md
    genReadme(projectOut, args.projName)
    
if __name__ == "__main__":
    args = parse_args()
    genCppProject(args)
