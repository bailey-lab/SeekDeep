#!/usr/bin/env python

import os, sys, glob, shutil
import argparse

def rename_all(debug = False):
    for i, fn in enumerate(sorted(glob.glob('*'))):
        ext = os.path.splitext(fn)[1]
        nfn = str(i) + ext
        print fn, '-->', nfn
        if not debug:
            shutil.move(fn, nfn)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action="store_true", default=False)
    return parser.parse_args()

def main():
    args = parse_args()
    rename_all(args.debug)

if __name__ == "__main__":
    sys.exit(main())
