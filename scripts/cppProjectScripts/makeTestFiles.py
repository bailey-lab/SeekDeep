#!/usr/bin/env python

import os, glob, sys, errno

def populateFile(fnp, testfnp):
    d = open(fnp).read().split('\n')
    classes = set()
    for line in d:
        if line.startswith("class "):
            t = line.rstrip().split(' ')
            classes.add(t[1])
    t = fnp.split('/src/')
    inc = '//#include "../src/' + t[1] + '"\n'
    with open(testfnp, 'w') as f:
        f.write('#include <catch.hpp>\n\n')
        f.write(inc)
        f.write('//using namespace bib;\n\n')
        for c in classes:
            f.write('//TEST_CASE("Basic tests for ')
            f.write(c + '", "['+c+']" ){}\n\n')
        if not classes:
            c = os.path.basename(fnp).replace(".hpp", "")
            f.write('//TEST_CASE("Basic tests for ')
            f.write(c + '", "['+c+']" ){}\n\n')

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def touch(fname, times=None):
    with file(fname, 'a'):
        os.utime(fname, times)

root_path = os.path.dirname(os.path.abspath(__file__))
root_path = os.path.join(root_path, "../")
src_folders = glob.glob("{d}/src".format(d=root_path))

incs = []

for path in src_folders:
    for root, dirs, files in os.walk(path):
        for fn in files:
            fnp = os.path.join(root, fn)
            fnp = os.path.abspath(fnp)
            if not fnp.endswith(".hpp"):
                continue
            testfnp = fnp.replace("/src/", "/test/src/").replace(".hpp", "Tester.cpp")
            d = os.path.dirname(testfnp)
            if not os.path.exists(d):
                mkdir_p(d)
            if not os.path.exists(testfnp):
                touch(testfnp)
                populateFile(fnp, testfnp)
