#!/usr/bin/python

import os, glob, sys, re

def process(fnp):
    fn, ext = os.path.splitext(fnp)
    fn_cpp = os.path.basename(fn) + ".cpp"

    s = open(fnp, 'r').read()
    match = re.search(r"""(?s)(.*
class.*\};)
(.*)""", s)
    if match:
        header = match.group(1)
        cpp = match.group(2)
    else:
        return

    with open(fnp, 'w') as f:
        f.write(header)
    with open(fn + ".cpp", 'w') as f:
        f.write('#include "' + os.path.basename(fnp) + '"\n')
        f.write(cpp)

#d = os.path.dirname(os.path.abspath(__file__))
#d = os.path.join(d, "../")
d = "/home/mjp//njh-cpp"
src_folders = glob.glob("{d}/src".format(d=d))

for path in src_folders:
    for root, dirs, files in os.walk(path):
        for fn in files:
            fnp = os.path.join(root, fn)
            fnp = os.path.abspath(fnp)
            if fnp.endswith(".h") or fnp.endswith(".hpp"):
                print(fnp)
                process(fnp)
