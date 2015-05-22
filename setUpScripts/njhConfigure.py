#!/usr/bin/env python

import shutil, os, argparse, sys, stat
from setUpScripts.utils import Utils
from setUpScripts.genFuncs import genHelper

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-prefix', type=str, nargs=1)
    parser.add_argument('-externalLibDir', type=str, nargs=1)
    parser.add_argument('-CC', type=str, nargs=1)
    parser.add_argument('-CXX', type=str, nargs=1)
    parser.add_argument('-libs', type=str)
    parser.add_argument('-name', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    prefix = "";
    external = "external";
    CC = genHelper.determineCC(args)
    CXX = genHelper.determineCXX(args)
    if(args.externalLibDir):
        external = args.externalLibDir[0];
    
    if(args.prefix):
        prefix = args.prefix[0];
        cmd = "setUpScripts/generateCompFile.py -prefix "  + prefix + " -installName " + args.name + " -outFilename compfile.mk -externalLoc " + external + " -CC " + CC + "  -CXX " + CXX + " -neededLibs " + args.libs + " -outname " + args.name + " "

    else:
        cmd = "setUpScripts/generateCompFile.py "  + " -installName " + args.name + " -outFilename compfile.mk -externalLoc " + external + " -CC " + CC + "  -CXX " + CXX  + "  -neededLibs " + args.libs + " -outname " + args.name + " "

    
    Utils.run(cmd)
    
main()
