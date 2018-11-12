#!/usr/bin/env python2

import shutil, os, argparse, sys, stat, time

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "setUpScripts"))
from genFuncs import genHelper
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyUtils"))
from utils import Utils
from color_text import ColorText as CT


def genArgParsePythonCompletes(programNames):
    ret = """
_argParsePys()
{
    local cur prev opts base
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    if [[ ${cur} == -* ]]; then
        opts=$(for x in `${COMP_WORDS[0]} -h | grep " -" | sed "s/^. *-/-/g" | sed "s/   .*//g" | sed "s/, / /g"`; do echo ${x} ; done )
        COMPREPLY=($(compgen -W "${opts}" -- ${cur}))
    else
        _filedir
    fi
   return 0
}
"""
    for name in programNames:
        ret += "complete -F _argParsePys {programName}\n".format(programName = name)
    return ret

def genSetUpPyCompletes():
    setUpPrograms = ["setup.py", "mapSrc.py", "needToRecompile.py", "fileModAffect.py", "configure.py"]
    return genArgParsePythonCompletes(setUpPrograms)


def addSetUpPyCompletes(dest, outFilename):
    with open(os.path.join(dest,"bash_completion.d",outFilename), "w") as f:
        f.write(genSetUpPyCompletes())

def genMultiRingBashCompleteStr(programNames):
    ret =  """
_njhCppTools()
{
    local cur prev opts base
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    if [[ $COMP_CWORD -lt 2 ]] ; then
        opts=$(for x in `${COMP_WORDS[0]} | grep ")" | sed "s/.*) //g"`; do echo ${x} ; done )
        COMPREPLY=($(compgen -W "${opts}" -- ${cur}))
    elif [[ ${cur} == -* ]]; then
        if [[ ${COMP_WORDS[1]} == batch* ]]; then
            rest="${COMP_WORDS[@]:1:${#COMP_WORDS[@]} }"
            if [[ $rest != *"-getFlags"* ]]; then
                  rest="$rest -getFlags"
            fi
            newopts=$(${COMP_WORDS[0]} $rest | column -t | cut -f 1 -d " " | cut -f 1 -d ,)
            COMPREPLY=( $(compgen -W "${newopts}" -- ${cur}) )
        else
            newopts=$(${COMP_WORDS[0]} ${COMP_WORDS[1]} -getFlags | column -t | cut -f 1 -d " " | cut -f 1 -d ,)
            COMPREPLY=( $(compgen -W "${newopts}" -- ${cur}) )
        fi
    else
        if [[ ${prev} == -run ]]; then
            opts=$(for x in `${COMP_WORDS[0]} | grep ")" | sed "s/.*) //g"`; do echo ${x} ; done )
            COMPREPLY=($(compgen -W "${opts}" -- ${cur}))
        else
            _filedir
        fi
    fi
   return 0
}
"""
    for name in programNames:
        ret += "complete -F _njhCppTools {programName}\n".format(programName = name)
    return ret

def addMultiRingComletes(dest,programNames, outFilename):
    with open(os.path.join(dest,"bash_completion.d",outFilename), "w") as f:
        f.write(genMultiRingBashCompleteStr(programNames))
        
def genSingleCmdBashCompleteStr(programNames):
    ret =  """
_singleNJHCppTools()
{
    local cur prev opts base
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    if [[ ${cur} == -* ]]; then
        newopts=$(${COMP_WORDS[0]} -getFlags | column -t | cut -f 1 -d " " | cut -f 1 -d ,)
        COMPREPLY=( $(compgen -W "${newopts}" -- ${cur}) )
    else
        _filedir
    fi
   return 0
}

"""
    for name in programNames:
        ret += "complete -F _singleNJHCppTools {programName}\n".format(programName = name)
    return ret

def addSingleCmdComletes(dest,programNames, outFilename):
    with open(os.path.join(dest,"bash_completion.d",outFilename), "w") as f:
        f.write(genSingleCmdBashCompleteStr(programNames))

def genBashCompleteFolder(dest):
    os.mkdir(os.path.join(dest, "bash_completion.d"))

def fileInfoHeader(headerName, author):
    
    return """
//  {headerName}
//
//  Created by {author} on {date}.
//  Copyright (c) {year} {author}. All rights reserved.
//
""".format(headerName=headerName, author=author,year=time.strftime("%Y"),date=time.strftime("%Y/%m/%d"))

def startHeader(headerName, author):
    
    return """#pragma once
//
""" + fileInfoHeader(headerName, author)

def startCpp(nameStub, author):
    return fileInfoHeader(nameStub + ".cpp", author) + """
    
#include "{name}.hpp"
    
    """.format(name = nameStub)



def genRing(runnerName,ringName, projNamespace, internalIncludes, externalIncludes, parentPath, dest, author, placeHolderFunc):
    externalIncludes = ["njhcpp.h"] + externalIncludes
    if not parentPath.endswith("/"):
        parentPath += "/"
    ringDestName = os.path.join(dest, ringName)
    if os.path.exists(ringDestName) or os.path.exists(ringDestName + ".h"):
        print "Error, " + ringDestName + " already exists"
        exit(1)
    #create main dir

    os.mkdir(ringDestName)
    #create main header to include the ring
    with open(ringDestName + ".h", "w") as f:
        mainHeaderOut = startHeader( ringName + ".h", author) + """
#include "{parentPath}{name}/{name}SetUp.hpp"
#include "{parentPath}{name}/{name}Runner.hpp"
""".format(name=ringName,parentPath = parentPath)
        f.write(mainHeaderOut)

    #create setUp header
    with open(os.path.join(ringDestName,ringName + "SetUp.hpp"), "w") as f:
        defaultHeader = startHeader(ringName + "SetUp.hpp", author)
        for eInclude in externalIncludes:
            defaultHeader += "#include <" + eInclude + ">\n"
        for iInclude in internalIncludes:
            defaultHeader += "#include \"" + iInclude + "\"\n"
        defaultHeader += """
namespace {projNamespace} {{

class {name}SetUp : public njh::progutils::programSetUp {{

 public:
    using programSetUp::programSetUp; //include programSetUp's constructors
}};
}} // namespace {projNamespace}
""".format(name =ringName, projNamespace = projNamespace)
        f.write(defaultHeader)
    #create setUp cpp
    with open(os.path.join(ringDestName,ringName + "SetUp.cpp"), "w") as f:
        infoHeader = startCpp(ringName + "SetUp", author)
        infoHeader +="""
namespace {projNamespace} {{

}} // namespace {projNamespace}
""".format(projNamespace = projNamespace)
        f.write(infoHeader)
    #create runner header
    with open(os.path.join(ringDestName,ringName + "Runner.hpp"), "w") as f:
        infoHeader = startHeader(ringName + "Runner.hpp", author)
        infoHeader +="""
#include "{name}SetUp.hpp"

namespace {projNamespace} {{

class {name}Runner : public njh::progutils::programRunner {{
 public:
  {name}Runner();
  
  static int {placeHolderFunc}(std::map<std::string, std::string> inputCommands);

}};
}} // namespace {projNamespace}
""".format(name = ringName,projNamespace = projNamespace, placeHolderFunc= placeHolderFunc)
        f.write(infoHeader)
    #create runner cpp
    with open(os.path.join(ringDestName,ringName + "Runner.cpp"), "w") as f:
        infoHeader = startCpp(ringName + "Runner", author)
        infoHeader +="""
namespace {projNamespace} {{

{name}Runner::{name}Runner()
    : njh::progutils::programRunner({{addFunc("{placeHolderFunc}", {placeHolderFunc}, false)}},
                    "{runnerName}") {{}}
                    
int {name}Runner::{placeHolderFunc}(std::map<std::string, std::string> inputCommands) {{
  {name}SetUp setUp(inputCommands);
  std::string name = "World";
  setUp.setOption(name, "--name", "Someone\'s Name");
  setUp.finishSetUp(std::cout);
  std::cout << "From {name} {placeHolderFunc}, Hello " << name << "!" << std::endl;
  return 0;
}}
                    
}} // namespace {projNamespace}
""".format(name = ringName, projNamespace = projNamespace, placeHolderFunc= placeHolderFunc, runnerName = runnerName)
        f.write(infoHeader)

def genOneRing(runnerName, ringName, projNamespace, rings, externalIncludes, parentPath, dest, author, placeHolderFunc):
    externalIncludes = ["njhcpp.h"] + externalIncludes
    if not parentPath.endswith("/"):
        parentPath += "/"
    ringDestName = os.path.join(dest, ringName)
    if os.path.exists(ringDestName) or os.path.exists(ringDestName + ".h"):
        print "Error, " + ringDestName + " already exists"
        exit(1)
    #create main dir

    os.mkdir(ringDestName)
    #create main header to include the ring
    with open(ringDestName + ".h", "w") as f:
        mainHeaderOut = startHeader( ringName + ".h", author) + """
#include "{parentPath}{name}/{name}SetUp.hpp"
#include "{parentPath}{name}/{name}Runner.hpp"
""".format(name=ringName,parentPath = parentPath)
        f.write(mainHeaderOut)

    #create setUp header
    with open(os.path.join(ringDestName,ringName + "SetUp.hpp"), "w") as f:
        defaultHeader = startHeader(ringName + "SetUp.hpp", author)
        for eInclude in externalIncludes:
            defaultHeader += "#include <" + eInclude + ">\n"
        defaultHeader += """
namespace {projNamespace} {{

class {name}SetUp : public njh::progutils::programSetUp {{

 public:
    using programSetUp::programSetUp; //include programSetUp's constructors
}};
}} // namespace {projNamespace}
""".format(name =ringName, projNamespace = projNamespace)
        f.write(defaultHeader)
    #create setUp cpp
    with open(os.path.join(ringDestName,ringName + "SetUp.cpp"), "w") as f:
        infoHeader = startCpp(ringName + "SetUp", author)
        infoHeader +="""
namespace {projNamespace} {{

}} // namespace {projNamespace}
""".format(projNamespace = projNamespace)
        f.write(infoHeader)
    #create runner header
    with open(os.path.join(ringDestName,ringName + "Runner.hpp"), "w") as f:
        infoHeader = startHeader(ringName + "Runner.hpp", author)
        infoHeader +="""
#include "{name}SetUp.hpp"

namespace {projNamespace} {{

class {name}Runner : public njh::progutils::oneRing {{
 public:
  {name}Runner();
  
  static int {placeHolderFunc}(std::map<std::string, std::string> inputCommands);

}};
}} // namespace {projNamespace}
""".format(name = ringName,projNamespace = projNamespace, placeHolderFunc = placeHolderFunc)
        f.write(infoHeader)
    #create runner cpp
    with open(os.path.join(ringDestName,ringName + "Runner.cpp"), "w") as f:
        infoHeader = startCpp(ringName + "Runner", author)
        infoHeader +="\n"
        for iInclude in rings:
            infoHeader += "#include \"{parentPath}{addRingPrefix}/{addRing}.hpp\"\n".format(parentPath = parentPath, addRingPrefix = iInclude.replace("Runner", ""), addRing =iInclude )
        infoHeader +="""
namespace {projNamespace} {{

{name}Runner::{name}Runner()
    : njh::progutils::oneRing({{"""

        for ring in rings:
            infoHeader += "addRing(std::make_shared<{ring}>()),".format(ring = ring)
        infoHeader += """}},{{addFunc("{placeHolderFunc}", {placeHolderFunc}, false)}},
                    "{runnerName}") {{}}
                    
int {name}Runner::{placeHolderFunc}(std::map<std::string, std::string> inputCommands) {{
  {name}SetUp setUp(inputCommands);
  std::string name = "World";
  setUp.setOption(name, "--name", "Someone\'s Name");
  setUp.finishSetUp(std::cout);
  std::cout << "From {name} {placeHolderFunc}, Hello " << name << "!" << std::endl;
  return 0;
}}
                    
}} // namespace {projNamespace}
"""
        infoHeader = infoHeader.format(name = ringName, projNamespace = projNamespace,placeHolderFunc= placeHolderFunc,  runnerName = runnerName)
        f.write(infoHeader)


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

def genProgramHeader(dest, projName):
    with open(os.path.join(dest, "src/" + projName + "/programs.h"), "w") as f:
        headers = [progFile for progFile in os.listdir(os.path.join(dest, "src/" + projName + "/programs")) if progFile.endswith(".h")]
        f.write("#pragma once\n")
        f.write("\n")
        f.write("// Created on "+ time.strftime("%Y/%m/%d") +"\n")
        f.write("// Including headers in common\n")
        f.write("\n")
        for head in headers:
            f.write("#include \""+ projName + "/programs/" + head + "\"\n")
        
def genWholeProjInclude(dest, projName, addProgramDir):
     with open(os.path.join(dest, "src/" + projName + ".h"), "w") as f:
         f.write("#pragma once\n")
         f.write("\n")
         f.write("// Created on "+ time.strftime("%Y/%m/%d") +"\n")
         f.write("// Including whole project\n")
         f.write("\n")
         f.write("#include \""+ projName + "/common.h\"\n")
         if(addProgramDir):
             f.write("#include \""+ projName + "/programs.h\"\n")
         
 
def genMain(dest, projName):
    with open(os.path.join(dest, "src/" + "main.cpp"), "w") as f:
        f.write("\n")
        f.write("// Created on "+ time.strftime("%Y/%m/%d") +"\n")
        f.write("// main.cpp\n")
        f.write("\n")
        f.write("#include \""+ projName + ".h\"")
        f.write("\n")
        f.write("int main(int argc, char* argv[]){\n")
        f.write("  " + projName + "::" +  projName +"ProgramRunner runner;\n")
        f.write("  return runner.run(argc, argv);\n")
        f.write("}\n")

def genSingCmdMain(dest, projName):
    with open(os.path.join(dest, "src/" + "main.cpp"), "w") as f:
        f.write("\n")
        f.write("// Created on "+ time.strftime("%Y/%m/%d") +"\n")
        f.write("// main.cpp\n")
        f.write("\n")
        f.write("#include \""+ projName + ".h\"\n")
        f.write("#include <njhcpp.h> \n")
        f.write("\n")
        f.write("int main(int argc, char* argv[]){\n")
        mainContent = """
    njh::progutils::programSetUp setUp(argc, argv);
    std::string name = "World";
    setUp.setOption(name, "--name", "Name to say hello to", false);
    setUp.finishSetUp(std::cout);
    std::cout << "Hello " << name << "!\\n";
    return 0;
"""
        f.write(mainContent)
        f.write("}\n")
        


def genSrcSingleRingProgram(dest, projName, includes, externalIncludes, author):
    os.mkdir(os.path.join(dest, "src"))
    os.mkdir(os.path.join(dest, "src/" + projName))
    os.mkdir(os.path.join(dest, "src/" + projName + "/common"))
    os.mkdir(os.path.join(dest, "src/" + projName + "/programs"))     
    genCommon(dest, projName,  includes)
    
    genWholeProjInclude(dest, projName, True)
    intIncForRing = [projName + "/common.h"]
    genRing(projName, projName + "Program", projName, intIncForRing, externalIncludes, projName + "/programs", projName + "/" + "src/" + projName + "/programs", author, "hellowWorld")
    genMain(dest, projName)
    genProgramHeader(dest, projName)
    genBashCompleteFolder(projName)
    addMultiRingComletes(projName, [projName], projName)
    addSetUpPyCompletes(projName, "pyCompletes")
    

def genSrcWithOneRingProgram(dest, projName, includes, externalIncludes, author):
    os.mkdir(os.path.join(dest, "src"))
    os.mkdir(os.path.join(dest, "src/" + projName))
    os.mkdir(os.path.join(dest, "src/" + projName + "/common"))
    os.mkdir(os.path.join(dest, "src/" + projName + "/programs"))     
    genCommon(dest, projName,  includes)
    genWholeProjInclude(dest, projName, True)
    intIncForRing = [projName + "/common.h"]
    genOneRing(projName, projName + "Program", projName, [projName + "Sub1" + "Runner", projName + "Sub2" + "Runner"], externalIncludes, projName + "/programs", projName + "/" + "src/" + projName + "/programs", author, "hellowWorldMain")
    genRing(projName + "Sub1", projName + "Sub1", projName, intIncForRing, externalIncludes, projName + "/programs", projName + "/" + "src/" + projName + "/programs", author, "hellowWorld1")
    genRing(projName + "Sub2", projName + "Sub2", projName, intIncForRing, externalIncludes, projName + "/programs", projName + "/" + "src/" + projName + "/programs", author, "hellowWorld2")
    genProgramHeader(dest, projName)
    genMain(dest, projName)
    genBashCompleteFolder(projName)
    addMultiRingComletes(projName, [projName], projName)
    addSetUpPyCompletes(projName, "pyCompletes")

def genSrcWithOneCmdProgram(dest, projName, includes, externalIncludes, author):
    os.mkdir(os.path.join(dest, "src"))
    os.mkdir(os.path.join(dest, "src/" + projName))
    os.mkdir(os.path.join(dest, "src/" + projName + "/common"))
    os.mkdir(os.path.join(dest, "src/" + projName + "/programs"))     
    genCommon(dest, projName,  includes)
    genWholeProjInclude(dest, projName, False)
    genSingCmdMain(dest, projName)
    genBashCompleteFolder(projName)
    addSingleCmdComletes(projName, [projName], projName)
    addSetUpPyCompletes(projName, "pyCompletes")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-projName', type=str, nargs=1, required=True)
    parser.add_argument('-dest', type=str, nargs=1, required=True)
    parser.add_argument('-CC', type=str, nargs=1)
    parser.add_argument('-CXX', type=str, nargs=1)
    parser.add_argument('-externalLoc', type=str, nargs=1)
    parser.add_argument('-neededLibs', type=str, nargs=1)
    parser.add_argument('-author', type=str, required=True)
    parser.add_argument('-programType', type=str, required=True)
    return parser.parse_args()

def main():
    args = parse_args()
    externalIncludes = []
    stdLibraryInc = ["iostream", "string", "unistd.h", "vector", "cstdint", "cstdio", "cstddef", "utility", "map", "unordered_map", "algorithm"]
    projectOut = os.path.join(args.dest[0], args.projName[0])
    os.mkdir(projectOut)
    if args.programType == "singleRing":
        genSrcSingleRingProgram(projectOut, args.projName[0], stdLibraryInc, externalIncludes, args.author)
    elif args.programType == "oneRing":
        genSrcWithOneRingProgram(projectOut, args.projName[0], stdLibraryInc, externalIncludes, args.author)
    elif args.programType == "oneCmd":
        genSrcWithOneCmdProgram(projectOut, args.projName[0], stdLibraryInc, externalIncludes, args.author)
    else:
        raise Exception("Error, only singleRing, oneRing,oneCmd available for options to programType, was given " + args.programType )
    
    CC = genHelper.determineCC(args)
    CXX = genHelper.determineCXX(args)
    external = "external"
    outname = args.projName[0]
    prefix = "./"
    installName = args.projName[0]
    neededLibs = ["njhcppdev"]        
    if args.externalLoc:
        external = os.path.realpath(args.externalLoc[0])
    if args.neededLibs:
        neededLibs = ["njhcppdev"] + args.neededLibs[0].split(",")
    genHelper.generateCompfileFull(os.path.join(projectOut, "compfile.mk"), external, CC, CXX, outname, installName, prefix, neededLibs)
    with open(os.path.join(projectOut, "configure.py"), "w") as configFile:
        if(args.neededLibs):
            configFile.write(genHelper.mkConfigFileStr(outname, ",".join(neededLibs)))
        else:
            configFile.write(genHelper.mkConfigFileStr(outname, "njhcppdev"))
    os.chmod(os.path.join(projectOut, "configure.py"), stat.S_IXGRP | stat.S_IXOTH | stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IWUSR)
    exFrom = os.path.abspath(os.path.dirname(__file__))
    cpSetUpCmd = exFrom + "/copySetUpFiles.py -from " + exFrom +"/../../ -to " + projectOut
    print CT.boldBlack(cpSetUpCmd)
    Utils.run(cpSetUpCmd)
    cpMakefilesCmd = "cp " + exFrom + "/../cppSetUpFiles/*akefile* " + projectOut
    print CT.boldBlack(cpMakefilesCmd)
    Utils.run(cpMakefilesCmd)
    
main()
