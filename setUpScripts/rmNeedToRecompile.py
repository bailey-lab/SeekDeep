#!/usr/bin/env python

import fnmatch, subprocess, sys, os, argparse, re, copy

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
    parser.add_argument('-obj', type=str, nargs=1, required =True);
    return parser.parse_args()
def main():
    args = parse_args()
    allFiles = fileCollection.getAllSourceFiles(args.src[0])
    objectFiles = fileCollection.getObjectFiles(args.obj[0])
    srcFileDict = {}
    graph = headInGraph();
    objectFilesTrimedDic = {}
    for obj in objectFiles:
        srcName = obj.find(args.src[0])
        objectFilesTrimedDic[obj[srcName:]] = obj
        srcFileDict[os.path.basename(obj).replace(".o", "_cpp")] = obj
    
    for file in allFiles:
       
        statbuf = os.stat(file)
        if(".h" not in file):
            graph.addNode(os.path.basename(file).replace(".", "_"), fileNode.cppColor, "internal", statbuf.st_mtime,os.path.getsize(file))
            lastPeriod = file.rfind(".c")
            objNameFile = file[:lastPeriod] + ".o"
            if objNameFile in objectFilesTrimedDic.keys():
                statbuf = os.stat(objectFilesTrimedDic[objNameFile])
                graph.addObjecTime(os.path.basename(file).replace(".", "_"), statbuf.st_mtime)
        else:
            graph.addNode(os.path.basename(file).replace(".", "_"), fileNode.headerColor, "internal", statbuf.st_mtime, os.path.getsize(file))
        
    pattern = re.compile("#include.*\".*\.h")
    for file in allFiles:
        for i, line in enumerate(open(file)):
            for match in re.finditer(pattern, line):
                if(".h" in file): 
                    graph.addPair(os.path.basename(file).replace(".", "_"), os.path.basename(re.findall('"([^"]*)"', line)[0]).replace(".", "_"), fileNode.headToHeadColor)
                else:
                    graph.addPair(os.path.basename(file).replace(".", "_"), os.path.basename(re.findall('"([^"]*)"', line)[0]).replace(".", "_"), fileNode.cppToHeaderColor)
    graph.reset()
    needToRecompile = [];
    for nPos in range(len(graph.nodes_)):
        if not (graph.nodes_[nPos].value_.endswith("_c") or graph.nodes_[nPos].value_.endswith("_cpp")):
            for child in graph.getChildrenList(nPos):
                if not (graph.nodes_[child].value_.endswith("_h") or graph.nodes_[child].value_.endswith("_hpp") or child == nPos):
                    if graph.nodes_[nPos].modTime_ > graph.nodes_[child].objectModTime_:
                        needToRecompile.append(graph.nodes_[child].value_)
                    else:
                        graph.nodes_[child].visited_ = False
                else:
                    graph.nodes_[child].visited_ = False
        
    for need in needToRecompile:
        if(need in srcFileDict.keys()):
            os.remove(srcFileDict[need]);

main()