#!/usr/bin/env python3

import fnmatch, subprocess, sys, os, argparse, re, copy
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyUtils"))
from headInGraph import fileCollection
from headInGraph import fileNode
from headInGraph import headInGraph





'''
#print "source files"
for f  in allFiles:
    statbuf = os.stat(f)
    #print f + " " + str(statbuf.st_mtime)
#print "Object files: " 
for f  in objectFiles:
    statbuf = os.stat(f)
    #print f + " " + str(statbuf.st_mtime)
    '''
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-src', type=str, nargs=1, required =True);
    parser.add_argument('-outMain', type=str, nargs=1, required = True);
    parser.add_argument('-addSystem', dest = 'addSystem', action = 'store_true' );
    return parser.parse_args()
def main():
    args = parse_args()
    allFiles = fileCollection.getAllSourceFiles(args.src[0])
    
    graph = headInGraph();
    sizes = []
    for file in allFiles:
        sizes.append(os.path.getsize(file))
   # for s in sizes:
        #print max(1, int(float(s)/max(sizes) * 100))
        
    for file in allFiles:
        #print os.path.getsize(file)
        statbuf = os.stat(file)
        #print int(float(os.path.getsize(file))/max(sizes) * 100)
        if(".h" not in file):
            graph.addNode(os.path.basename(file).replace(".", "_"), fileNode.cppColor, "internal", statbuf.st_mtime, max(1, int(float(os.path.getsize(file))/max(sizes) * 50) ))
        else:
            graph.addNode(os.path.basename(file).replace(".", "_"), fileNode.headerColor, "internal", statbuf.st_mtime, max(1, int(float(os.path.getsize(file))/max(sizes) * 50)) )
    pattern = re.compile("^[\w]*#include.*\".*\.h")
    patternSystem = re.compile("^[\w]*#include.*\<.*\>")
    for file in allFiles:
        for i, line in enumerate(open(file)):
            for match in re.finditer(pattern, line):
                if(".h" in file): 
                    graph.addPair(os.path.basename(file).replace(".", "_"), os.path.basename(re.findall('"([^"]*)"', line)[0]).replace(".", "_"), fileNode.headToHeadColor)
                else:
                    graph.addPair(os.path.basename(file).replace(".", "_"), os.path.basename(re.findall('"([^"]*)"', line)[0]).replace(".", "_"), fileNode.cppToHeaderColor)
            if args.addSystem:
                for match in re.finditer(patternSystem, line):
                    modifiedHeaderName = ((line[(line.find("<") + 1):line.find(">")]).replace(".", "_")).replace("/", "__");
                    if modifiedHeaderName not in list(graph.nodePositions_.keys()):
                        graph.addNode(modifiedHeaderName, fileNode.externalHeaderColor, "external", 0, 1)
                    if(".h" in file): 
                        graph.addPair(os.path.basename(file).replace(".", "_"), modifiedHeaderName, fileNode.headToHeadColor)
                    else:
                        graph.addPair(os.path.basename(file).replace(".", "_"), modifiedHeaderName, fileNode.cppToHeaderColor)
    
    outMainFile = open(args.outMain[0], "w");
    graph.printGraphViz(outMainFile, "all", args.addSystem)

main()