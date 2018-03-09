#!/usr/bin/env python3


import shutil, os, argparse, sys
from shutil import ignore_patterns



def copyDir(src, dist, overWrite):
    if(os.path.isdir(dist)):
        if(overWrite):
            shutil.rmtree(dist)
        else:
            print("Error, directory " + str(dist) + " already exist, set overwrite to overwrite dir")
            exit(1)
    shutil.copytree(src, dist, ignore=ignore_patterns('*.pyc', '*~', '.*'))


def copySetUp(cwDir, distDir, overWrite=False):
    if(os.path.isfile(os.path.join(distDir, "setup.py"))):
        if(overWrite):
            shutil.copy(os.path.join(cwDir, "setup.py"), os.path.join(distDir, "setup.py"))
        else:
            print("Error, file " + os.path.join(distDir, "setup.py") + " already exist, set overwrite to overwrite setup.py")
            exit(1)
    else:
        shutil.copy(os.path.join(cwDir, "setup.py"), os.path.join(distDir, "setup.py"))
    copyDir(os.path.join(cwDir + "scripts"), os.path.join(distDir,"scripts"), overWrite)
    

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-from',dest = "fromDir", type=str, nargs=1, required =True)
    parser.add_argument('-to', type=str, nargs=1, required =True)
    parser.add_argument('-overWrite', dest = 'overWrite', action = 'store_true' )
    return parser.parse_args()

def main():
    args = parse_args()
    for t in args.to[0].split(","):
        copySetUp(args.fromDir[0],t, args.overWrite)
    
main()