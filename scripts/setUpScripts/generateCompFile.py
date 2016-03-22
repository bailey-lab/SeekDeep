#!/usr/bin/env python

import shutil, os, argparse, sys, stat
from genFuncs import genHelper


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-CC', type=str, nargs=1)
    parser.add_argument('-CXX', type=str, nargs=1)
    parser.add_argument('-outname', type=str, nargs=1)
    parser.add_argument('-outFilename', type=str, nargs=1, required = True)
    parser.add_argument('-externalLoc', type=str, nargs=1)
    parser.add_argument('-prefix', type=str, nargs=1)
    parser.add_argument('-installName', type=str, nargs=1)
    parser.add_argument('-neededLibs', type=str, nargs=1)
    parser.add_argument('-ldFlags', type=str)
    return parser.parse_args()

def main():
    """@todo: Also add on adding CXXFLAGS or LDFLAGS and the such or look for environment ones"""
    args = parse_args()
    CC = genHelper.determineCC(args)
    CXX = genHelper.determineCXX(args)
    external = "external"
    outname = "out"
    prefix = "./"
    installName = "out"
    neededLibs = "none"
    ldFlags = ""
    if args.ldFlags and "" != args.ldFlags:
        ldFlags = args.ldFlags        
    if args.externalLoc:
        external = args.externalLoc[0]
    if args.outname:
        outname = args.outname[0]
    if args.installName:
        installName = args.installName[0]
    if args.prefix:
        prefix = args.prefix[0]
    if args.neededLibs:
        neededLibs = args.neededLibs[0].split(",")
    genHelper.generateCompfileFull(args.outFilename[0], external, CC, CXX, outname, installName, prefix, neededLibs, ldFlags)
    
main()
