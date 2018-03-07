#!/usr/bin/env python3

import fnmatch, subprocess, sys, os, argparse, re



class fileNode():
    headerColor = "#0571b0";
    headToHeadColor = "#92c5de";
    externalHeaderColor = "#af8dc3";
    cppColor = "#ca0020";
    cppToHeaderColor = "#f4a582";
    modColor = "#d01c8b";
    modEdgeColor = "#f1b6da";
    unModColor = "#4dac26";
    unModEdgeColor = "#b8e186";
    
    def __init__(self, name, color, type, modTime, filesize):
        self.value_ = name
        self.modTime_ = modTime
        self.objectModTime_ = 0
        self.childrenEdges_ = []
        self.color_ = color
        self.visited_ = False
        self.type_ = type
        self.filesize_ = filesize
    def changeObjectModTime(self, objModTime):
        self.objectModTime_ = objModTime

class fileEdge():
    
    def __init__(self, childPos, edgeColor):
        self.childPos_ = childPos;
        self.edgeColor_ = edgeColor;
        
    

class headInGraph():
    def __init__(self):
        self.nodes_ = []
        self.nodePositions_ = {}
    
    def addNode(self, name, color,type, modTime, filesize):
        self.nodePositions_[name] = len(self.nodes_)
        self.nodes_.append(fileNode(name,color,type, modTime, filesize))
    
    def addPair(self, including, beingIncluded, edgeColor):
        if beingIncluded not in self.nodePositions_:
            raise Exception("Header " + beingIncluded + " not found in source folder, being included by " + including)
        self.nodes_[self.nodePositions_[beingIncluded]].childrenEdges_.append(fileEdge(self.nodePositions_[including], edgeColor))
    
    def setNodeColorAll(self, newColor):
        for e, node in enumerate(self.nodes_):
            self.nodes_[e].color_ = newColor
    
    def setEdgeColorAll(self, newColor):
        for e, node in enumerate(self.nodes_):
            for child, node in enumerate(self.nodes_[e].childrenEdges_):
                self.nodes_[e].childrenEdges_[child].edgeColor_ = newColor
    
    def modChildren(self, nodePos, modColor, modEdgeColor):
        if not self.nodes_[nodePos].visited_:
            self.nodes_[nodePos].visited_ = True;
            self.nodes_[nodePos].color_ = modColor;
            for childEdgePos, child in enumerate(self.nodes_[nodePos].childrenEdges_) :
                self.nodes_[nodePos].childrenEdges_[childEdgePos].edgeColor_ = modEdgeColor
                self.modChildren(child.childPos_, modColor, modEdgeColor);
                
    def getChildrenList(self, nodePos):
        allChildren = []
        if not self.nodes_[nodePos].visited_:
            self.nodes_[nodePos].visited_ = True;
            allChildren.append(nodePos)
            for child in self.nodes_[nodePos].childrenEdges_:
                allChildren = allChildren + self.getChildrenList(child.childPos_)
        return allChildren
    
    def printInfo(self):
        for nPos in range(len(self.nodes_)):
            print(self.nodes_[nPos].value_)
            print(self.nodes_[nPos].visted_)
            if self.nodes_[nPos].visted_:
                print("\033[1;32mvisited true\033[0m")
            else:
                print("\033[1;31m visted is false\033[0m")
    
    def printChildren(self, nodePos, out):
        if not self.nodes_[nodePos].visited_:
            self.nodes_[nodePos].visited_ = True;
            out.write(self.nodes_[nodePos].value_ + " ");
            for child in self.nodes_[nodePos].childrenEdges_:
                self.printChildren(child.childPos_, out)
    
    def addObjecTime(self, name, objModTime):
        self.nodes_[self.nodePositions_[name]].objectModTime_ = objModTime
        
    def reset(self):
        for nPos in range(len(self.nodes_)):
            self.nodes_[nPos].visited_ = False
    
    def printGraphViz(self, out, title, outPutExternal, useFileSize =False):
        size = "1"
        if useFileSize:
            size = str( self.nodes_[nPos].filesize_) 
        out.write("digraph G  { \n")
        out.write("bgcolor =\"#000000\" \n")
        out.write("labelloc=\"t\" \n")
        out.write("fontcolor = \"#ffffff\"\n")
        out.write("fontsize = 20 \n")
        out.write("label = \"" + title + "\" \n")
        out.write("fixedsize = true; \n")
        for nPos in range(len(self.nodes_)):
            if(not outPutExternal and self.nodes_[nPos].type_ == "external"):
                continue
            out.write(self.nodes_[nPos].value_ + "[shape=circle,style=filled,fixedsize =false, color = \"#000000\", fillcolor =\"" + self.nodes_[nPos].color_ + "\", width = " + size + "]\n") 
        for nPos in range(len(self.nodes_)):
            if(not outPutExternal and self.nodes_[nPos].type_ == "external"):
                continue
            for child in self.nodes_[nPos].childrenEdges_:
                out.write(self.nodes_[nPos].value_ + " -> " + self.nodes_[child.childPos_].value_ + "[penwidth=5, color=\"") 
                out.write(child.edgeColor_);
                out.write("\"]\n")
        out.write("}\n")  
    
class fileCollection:
    @staticmethod
    def getCppFiles(searchDir):
        return [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(searchDir)
        for f in fnmatch.filter(files, '*.c*')]
    @staticmethod
    def getHeaderFiles(searchDir):
        return [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(searchDir)
        for f in fnmatch.filter(files, '*.h*')]
    @staticmethod
    def getObjectFiles(searchDir):
        return [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(searchDir)
        for f in fnmatch.filter(files, '*.o')]
    @staticmethod
    def getAllSourceFiles(searchDir):
        return fileCollection.getCppFiles(searchDir) + fileCollection.getHeaderFiles(searchDir);