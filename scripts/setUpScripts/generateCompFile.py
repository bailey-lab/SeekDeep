#!/usr/bin/env python2

import shutil, os, argparse, sys, stat
from genFuncs import genHelper


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-CC', type=str)
    parser.add_argument('-CXX', type=str)
    parser.add_argument('-outname', type=str)
    parser.add_argument('-outFilename', type=str, required = True)
    parser.add_argument('-externalLoc', type=str)
    parser.add_argument('-prefix', type=str)
    parser.add_argument('-installName', type=str)
    parser.add_argument('-neededLibs', type=str)
    parser.add_argument('-ldFlags', type=str)
    parser.add_argument('-cxxFlags', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    CC = genHelper.determineCC(args)
    CXX = genHelper.determineCXX(args)
    external = "external"
    outname = "out"
    prefix = "./"
    installName = "out"
    neededLibs = "none"
    ldFlags = ""
    cxxFlags = ""
    if args.ldFlags and "" != args.ldFlags:
        ldFlags = args.ldFlags  
    if args.cxxFlags and "" != args.cxxFlags:
        cxxFlags = args.cxxFlags        
    if args.externalLoc:
        external = args.externalLoc
    if args.outname:
        outname = args.outname
    if args.installName:
        installName = args.installName
    if args.prefix:
        prefix = args.prefix
    if args.neededLibs:
        neededLibs = args.neededLibs.split(",")
    genHelper.generateCompfileFull(args.outFilename, external, CC, CXX, outname, installName, prefix, neededLibs, ldFlags, cxxFlags)
    
main()
