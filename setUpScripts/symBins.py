#!/usr/bin/env python

import shutil, os, argparse, sys, stat, fnmatch

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-dest', type=str, nargs=1, required = True)
    parser.add_argument('-overWrite', dest = 'overWrite', action = 'store_true' )
    return parser.parse_args()

def main():
    args = parse_args()
    currentDirs = os.listdir(args.dest[0])
    for d in currentDirs:
        if os.path.isdir(d):
            nextDirs = fnmatch.filter(os.listdir(d), 'bin') 
            if nextDirs:
                for e in os.listdir(os.path.abspath(os.path.join(d,nextDirs[0]) )):
                    print e
                    print os.path.abspath(os.path.join(d,nextDirs[0], e) )
                    print os.path.abspath(os.path.join(args.dest[0], "bin", e) )
                    if(args.overWrite and os.path.lexists(os.path.abspath(os.path.join(args.dest[0], "bin", e) ))):
                        os.unlink(os.path.abspath(os.path.join(args.dest[0], "bin", e) ))
                    elif os.path.lexists(os.path.abspath(os.path.join(args.dest[0], "bin", e) )):
                        print os.path.abspath(os.path.join(args.dest[0], "bin", e) ) + " already exists, run with -overWrite to over write it"
                        exit(1)
                    os.symlink(os.path.abspath(os.path.join(d,nextDirs[0], e) ), os.path.abspath(os.path.join(args.dest[0], "bin", e) )) 
                        
    
main()
