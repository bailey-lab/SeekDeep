#!/usr/bin/env python2

import os, glob, sys

po = "#pragma once"

def is_header(fn):
    fn = fn.lower()
    if fn.endswith(".h") or fn.endswith(".hpp"):
        return True
    return False

def fixFile(fnp):
    q = open(fnp).read()
    with open(fnp, "w") as f:
        f.write(po + "\n")
        for line in q.split("\n"):
            if line.startswith("#define") or line.startswith("#ifndef"):
                print line
                continue
            if line.startswith("#endif"):
                print line
                continue
            f.write(line + "\n")

d = os.path.dirname(os.path.abspath(__file__))
d = os.path.join(d, "../")
src_folders = glob.glob("{d}/src".format(d=d))

for path in src_folders:
    for root, dirs, files in os.walk(path):
        for fn in files:
            if not is_header(fn):
                continue
            fnp = os.path.join(root, fn)
            fnp = os.path.abspath(fnp)
            firstline = open(fnp).readline().rstrip()
            if po != firstline:
                print "fixing", fnp
                fixFile(fnp)

