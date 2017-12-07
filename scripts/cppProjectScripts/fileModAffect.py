#!/usr/bin/env python2

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
    parser.add_argument('-header', type=str, nargs=1, required =True);
    parser.add_argument('-outMod', type=str, nargs=1, required = True);
    return parser.parse_args()
def main():
    args = parse_args()
    allFiles = fileCollection.getAllSourceFiles(args.src[0])
    
    graph = headInGraph();
    allFileSize = []
    for file in allFiles:
        allFileSize.append(os.path.getsize(file))
        statbuf = os.stat(file)
        if(".h" not in file):
            graph.addNode(os.path.basename(file).replace(".", "_"), fileNode.cppColor, "internal", statbuf.st_mtime, max(1, 50 * (os.path.getsize(file) /float(max(allFileSize)) ) ) )
        else:
            graph.addNode(os.path.basename(file).replace(".", "_"), fileNode.headerColor, "internal", statbuf.st_mtime,  max(1, 50 * (os.path.getsize(file) /float(max(allFileSize)) ) ) )
    pattern = re.compile("^[\w]*#include.*\".*\.h")
    for file in allFiles:
        for i, line in enumerate(open(file)):
            for match in re.finditer(pattern, line):
                if(".h" in file): 
                    graph.addPair(os.path.basename(file).replace(".", "_"), os.path.basename(re.findall('"([^"]*)"', line)[0]).replace(".", "_"), fileNode.headToHeadColor)
                else:
                    graph.addPair(os.path.basename(file).replace(".", "_"), os.path.basename(re.findall('"([^"]*)"', line)[0]).replace(".", "_"), fileNode.cppToHeaderColor)
    graph.setNodeColorAll(fileNode.unModColor)
    graph.setEdgeColorAll(fileNode.unModEdgeColor)
    outFileTest = open(args.outMod[0], "w")
    graph.modChildren(graph.nodePositions_[args.header[0]], fileNode.modColor, fileNode.modEdgeColor)
    graph.nodes_[graph.nodePositions_[args.header[0]]].color_ = "#ff0000"
    graph.printGraphViz(outFileTest, args.header[0], False)

main()