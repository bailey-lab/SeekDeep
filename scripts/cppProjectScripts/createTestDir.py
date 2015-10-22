#!/usr/bin/env python

import shutil, os, argparse, sys, stat,errno
import CppHeaderParser
from string import replace
from argparse import Action
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyUtils"))
from color_text import ColorText as CT
from headInGraph import *

testerBodyTemplate = """
/*
TEST_CASE("Basic tests for {REPLACETHIS}", "[{REPLACETHIS_DETAILED}]" ){{
  SECTION("GIVE SECTION NAME"){{
      YOUR CODE GOES HERE
        NORMALLY END WITH A REQUIRE STATEMENT e.g.
        REQUIRE(TESTVAL1 == YOURVAL);
  }}
}}
*/
"""

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
        
def mkdir_p_forFile(path):
    mkdir_p(os.path.dirname(path))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--src', type=str, required = True)
    parser.add_argument('--outDir', type=str, required = True)
    parser.add_argument("--overWrite", action = 'store_true')
    #parser.add_argument("--update", action = 'store_true')
    return parser.parse_args()

def getFuncDetailed(func):
    ret = ""
    ret = ret + (func["rtnType"] + " ")
    ret = ret + (func["name"] + " (")
    count = 0
    for par in func["parameters"]:
        if(count != 0):
            ret = ret + (",")
        count +=1
        ret = ret + (par["raw_type"])
        if(par["reference"]):
            ret = ret + ("&")
        elif par["pointer"]:
            ret = ret + ("*")
        ret = ret + (" " + par["name"])
    ret = ret + ")"
    return ret

def createTestMain(path, overWrite):
    mainBody = """
// based off https://github.com/philsquared/Catch/blob/master/docs/tutorial.md

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include <catch.hpp>

    """
    mainPath = os.path.join(path, "main.cpp")
    if os.path.exists(mainPath):
        if overWrite:
            os.remove(mainPath)
        else:
            print mainPath, "already exists, use --overWrite to remove current"
            return
    with open(mainPath, "w") as mainFile:
        mainFile.write(mainBody)
def copyMakefile(fromLoc, dest, overWrite):
    if os.path.exists(dest):
        if overWrite:
            os.remove(dest)
        else:
            print dest, "already exists, use --overWrite to replace it"
            return
    shutil.copy(fromLoc, dest)

def main():
    args = parse_args()
    headers = fileCollection.getHeaderFiles(args.src)
    
    for head in headers:
        try:
            cppHeader = CppHeaderParser.CppHeader(head)
        except CppHeaderParser.CppParseError as e:
            print(e)
            sys.exit(1)
            print CT.boldBlack("Class public methods")

        if(len(cppHeader.classes) + len(cppHeader.functions) > 0):
            testerCppPath = os.path.join(args.outDir,head.replace(".hpp", "Tester.cpp"))
            mkdir_p_forFile(testerCppPath)
            if os.path.exists(testerCppPath):
                if args.overWrite:
                    os.remove(testerCppPath)
                else:
                    print "Skipping", testerCppPath, "it already exist, use --overWrite to replace"
                    continue
            with open(testerCppPath, "w") as testerFile:
                testerFile.write("#include <catch.hpp>\n")
                testerFile.write("#include \"" + "../" + head + "\"\n")
                for func in cppHeader.functions:
                    testerFile.write(testerBodyTemplate.format(REPLACETHIS=func["name"], REPLACETHIS_DETAILED = getFuncDetailed(func)))
                for k in cppHeader.classes.keys():
                    for i in range(len(cppHeader.classes[k]["methods"]["public"])):
                        testerFile.write(testerBodyTemplate.format(REPLACETHIS=cppHeader.classes[k]["methods"]["public"][i]["name"], REPLACETHIS_DETAILED = getFuncDetailed(cppHeader.classes[k]["methods"]["public"][i])))
    createTestMain(os.path.join(args.outDir, args.src), args.overWrite)
    copyMakefile("scripts/cppSetUpFiles/unitTest/Makefile", os.path.join(args.outDir, "Makefile"), args.overWrite)
    return 0

main()
