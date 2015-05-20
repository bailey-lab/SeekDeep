#!/usr/bin/python

import os, glob, sys

def fixHeader(fnp,namepsace):
    q = open(fnp).read()
    with open(fnp, "w") as f:
        dump = False
        inc = False
        opened = False
        closed = False
        for line in q.split("\n"):
            if dump and inc:
                f.write(line + "\n")
                continue
            if dump:
                if line.startswith("#"):
                    f.write("\n} // namespace " + namespace + "\n\n")
                    closed = True
                    f.write(line + "\n")
                    inc = True
                    continue
                f.write(line + "\n")
                continue
            if line.startswith("#") or line.startswith("//") or line.strip() == "":
                f.write(line + "\n")
                continue
            dump = True
            f.write("namespace " + namespace + " {\n\n")
            opened = True
            f.write(line + "\n")
        if opened and not closed:
            f.write("\n} // namespace " + namespace + "\n")

def fixCPP(fnp, namepsace):
    q = open(fnp).read()
    with open(fnp, "w") as f:
        dump = False
        opened = False
        for line in q.split("\n"):
            if dump:
                f.write(line + "\n")
                continue
            if line.startswith("#") or line.startswith("/") or line.strip() == "":
                f.write(line + "\n")
                continue
            dump = True
            f.write("namespace " + namespace + " {\n\n")
            opened = True
            f.write(line + "\n")
        if opened:
            f.write("\n} // namespace " + namespace + "\n")

#d = os.path.dirname(os.path.abspath(__file__))
#d = os.path.join(d, "../")
#d = "/home/mjp/bib-cpp/"
#src_folders = glob.glob("{d}/src".format(d=d))

def checkN(fnp):
    s = open(fnp, "r").read()
    c = s.count("namespace")
    if 2 != c and 0 != c and 3 != c:
        print fnp, "bad count", c
        raise Exception(fnp)
    
def namepsaceSrcTree(sourceFolder):
    for root, dirs, files in os.walk(sourceFolder):
        for fn in files:
            if fn == "main.cpp":
                continue
            fnp = os.path.join(root, fn)
            fnp = os.path.abspath(fnp)
            if fnp.endswith(".h") or fnp.endswith(".hpp"):
                fixHeader(fnp)
            elif fnp.endswith(".cpp"):
                fixCPP(fnp)
            else:
                continue
            checkN(fnp)
