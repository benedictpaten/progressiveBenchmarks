#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Cactus progressive renames sequences to make sure they are compliant 
with mafjoin.  This is hopefully a temporary hack. This class goes through
a project and constructs a map of all the original sequence names to 
the output sequence names

"""

import os
import sys
import math
import xml.etree.ElementTree as ET
from optparse import OptionParser
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.preprocessor.cactus_addFastaHeaderDots import fixHeader
from cactus.preprocessor.cactus_preprocessor import fileList
from sonLib.bioio import fastaRead

class NamingMap:
    def __init__(self):
        self.nameMap = None
    
    def cactusName(self, sequenceName):
        assert sequenceName in self.nameMap
        return self.nameMap[sequenceName]
    
    def readProject(self, projectXmlPath):
        mcProj = MultiCactusProject()
        mcProj.readXML(projectXmlPath)
        mcProj.mcTree.nameUnlabeledInternalNodes()
        self.nameMap = dict()
        for leaf in mcProj.mcTree.getLeaves():
            eventName = mcProj.mcTree.getName(leaf)
            sequencePath = mcProj.sequencePath(eventName)
            for sequenceFile in fileList(sequencePath):
                if not os.path.isdir(sequenceFile):
                    self.processSequence(eventName, sequenceFile)
    
    def readExperiment(self, experimentXmlPath):
        expXml = ET.parse(experimentXmlPath).getroot()
        exp = ExperimentWrapper(expXml)
        self.nameMap = dict()
        for eventName, sequencePath in exp.seqMap.items():
            for sequenceFile in fileList(sequencePath):
                if not os.path.isdir(sequenceFile):
                    self.processSequence(eventName, sequenceFile)
    
    def processSequence(self, eventName, sequencePath):
        fileHandle = open(sequencePath, "r")
        for header, sequence in fastaRead(fileHandle):
            fixedHeader = fixHeader(header, event=eventName)
            if header in self.nameMap:
                assert self.nameMap[header] == fixedHeader
            else:
                self.nameMap[header] = fixedHeader

def applyNamingToMaf(experimentPath, inputPath, outputPath):
    inFile = open(inputPath, "r")
    outFile = open(outputPath, "w")
    nameMap = NamingMap()
    nameMap.readExperiment(experimentPath)
    
    for line in inFile:
        tokens = line.split()
        if len(tokens) >= 2 and tokens[0] == 's':
            name = tokens[1]
            outFile.write(line.replace(name, nameMap.cactusName(name), 1))
        else:
            outFile.write(line)
    
    inFile.close()
    outFile.close()
    
def main():
    usage = "usage: %prog <experiment> <input maf> <output maf>"
    description = "Apply naming conventions to maf"
    parser = OptionParser(usage=usage, description=description)
    
    options, args = parser.parse_args()
    
    if len(args) != 3:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    applyNamingToMaf(args[0], args[1], args[2])
    
if __name__ == '__main__':    
    main()