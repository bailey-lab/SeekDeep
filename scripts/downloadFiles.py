#!/usr/bin/env python

import shutil, os, argparse, sys, stat,urllib

def parse_args():
    parser = argparse.ArgumentParser(description="Take in a file where the first column holds the url of a file to be downloaded, will overwrite current files if they exist")
    parser.add_argument('-file', type=str, required = True)
    
    return parser.parse_args()

def main():
    args = parse_args()
    with open(args.file, "r") as f:
        for line in f:
            lineSplit = line.split()
            urllib.urlretrieve(lineSplit[0], os.path.basename(lineSplit[0]))
    
main()