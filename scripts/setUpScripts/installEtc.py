#!/usr/bin/env python

import shutil, os, argparse, sys, stat, fnmatch

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-etcFolder', type=str, required =True)
    parser.add_argument('-dest', type=str, required =True)
    parser.add_argument('-rmDir', action = "store_true", help = "Remove directory first")
    return parser.parse_args()

def main():
    args = parse_args()
    newDirName = os.path.join(args.dest,"etc")
    if os.path.exists(newDirName) and args.rmDir:
        shutil.rmtree(newDirName)
    elif os.path.exists(newDirName):
        raise Exception("Destination directory already exists, use -rmDir to overWrite")
    shutil.copytree(args.etcFolder, newDirName)
    
main()
