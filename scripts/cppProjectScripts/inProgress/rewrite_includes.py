#!/usr/bin/python

import os, glob, sys

def fixFile(fnp):
    q = open(fnp).read()
    with open(fnp, "w") as f:
        for line in q.split("\n"):
            if line.startswith("#include <seq") or line.startswith("#include <prog"):
                line = line.replace("<", '"').replace(">",'"')
            f.write(line + "\n")

#d = os.path.dirname(os.path.abspath(__file__))
#d = os.path.join(d, "../")
d = "/home/mjp/njh-cpp/"
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
            else:
                continue
