#!/usr/bin/env python3

import shutil, os, argparse, sys, stat, fnmatch

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-etcFolder', type=str, required =True)
    parser.add_argument('-dest', type=str, required =True)
    parser.add_argument('-rmDir', action = "store_true", help = "Remove directory first")
    return parser.parse_args()

def main():
    args = parse_args()
    allDirs = [x for x in os.walk(args.etcFolder[0])]
    installEtcDest = os.path.join(args.dest,"etc");
    for dInfo in allDirs:
        newDirName = dInfo[0].replace(args.etcFolder[0].rstrip('/'), args.installEtcDest.rstrip('/'))
        #print newDirName
        if not os.path.exists(newDirName):
            os.mkdir(newDirName)
        elif args.rmDir:
            shutil.rmtree(newDirName)
            os.mkdir(newDirName)
        for f in dInfo[2]:
            shutil.copy(os.path.join(dInfo[0], f), os.path.join(newDirName, f))
main()
