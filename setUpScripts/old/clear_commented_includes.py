#!/usr/bin/python

import os, glob, sys

def fixFile(fnp):
    q = open(fnp).read()
    with open(fnp, "w") as f:
        for line in q.split("\n"):
            if line.startswith("//#include"):
                continue
            f.write(line + "\n")

    q = open(fnp).read()
    with open(fnp, "w") as f:
        f.write(q.replace("\n\n\n", "\n"))

#d = os.path.dirname(os.path.abspath(__file__))
#d = os.path.join(d, "../")
d = "/home/mjp/bib-cpp/"
src_folders = glob.glob("{d}/src".format(d=d))

for path in src_folders:
    for root, dirs, files in os.walk(path):
        for fn in files:
            if fn == "main.cpp":
                continue
            fnp = os.path.join(root, fn)
            fnp = os.path.abspath(fnp)
            if fnp.endswith(".h") or fnp.endswith(".hpp") or fnp.endswith(".cpp"):
                fixFile(fnp)
