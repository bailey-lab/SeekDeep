#!/usr/bin/env python

import shutil, os, argparse, sys, stat
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyUtils"))

from utils import Utils
from genFuncs import genHelper

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-prefix', type=str)
    parser.add_argument('-externalLibDir', type=str)
    parser.add_argument('-CC', type=str, nargs = 1)
    parser.add_argument('-CXX', type=str, nargs = 1)
    parser.add_argument('-libs', type=str)
    parser.add_argument('-name', type=str, required = True)
    return parser.parse_args()

def main():
    args = parse_args()
    prefix = "";
    external = "external";
    CC = genHelper.determineCC(args)
    CXX = genHelper.determineCXX(args)
    if(args.externalLibDir):
        external = args.externalLibDir;
    cmd = """scripts/setUpScripts/generateCompFile.py -installName {name}
     -outFilename compfile.mk -externalLoc {external} -CC {CC} -CXX {CXX}
      -neededLibs {libs} -outname {name}"""
    if args.prefix and args.prefix != "":
        prefix = args.prefix;
        cmd += " -prefix {prefix}"
    cmd = " ".join(cmd.split())
    cmd = cmd.format(name = args.name, external = external, CC=CC, CXX=CXX, libs = args.libs, prefix = prefix)
    Utils.run(cmd)
    
main()
