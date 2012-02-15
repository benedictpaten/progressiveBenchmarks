#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" write a summary table in csv format for jobtree
and maf comparison results. 

"""

import os
import xml.etree.ElementTree as ET
import sys
import copy

from optparse import OptionParser
from progressiveBenchmarks.src.params import Params
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.experimentWrapper import ExperimentWrapper

class Summary:
    Header = Params.Header + \
    ["Run_Time", "Clock_Time", 
     "Base_tot_run", "Base_tot_clock", "Base_max_clock", 
     "Core1_tot_run", "Core1_tot_clock", "Core1_max_clock", 
     "Core2_tot_run", "Core2_tot_clock", "Core2_max_clock",
     "SelfBlast_tot_run", "SelfBlast_tot_clock", "SelfBlast_max_clock",
     "Blast_tot_run", "Blast_tot_clock", "Blast_max_clock",
     "Sensitivity", "Specificity", "Bal. Accuracy", 
     "Root Growth", "Avg Growth", "BBL_Min", "BBL_Max", "BBL_Avg", "TGS_Min", 
     "TGS_Max", "TGS_Avg"]
    SensIdx = len(Header)
    SpecIdx = SensIdx + 1 
    def __init__(self):
        self.table = []

    def addRow(self, catName, params, jobTreeStatsPath, mafCompPath, 
               treeStatsPath, projPath):
        if os.path.isfile(jobTreeStatsPath) and os.path.isfile(mafCompPath):
            mafXmlRoot = ET.parse(mafCompPath).getroot()
            jtXmlRoot = ET.parse(jobTreeStatsPath).getroot()
            tsXmlRoot = None
            if os.path.isfile(treeStatsPath):
                tsXmlRoot = ET.parse(treeStatsPath).getroot()
            if projPath is not None:
                project = MultiCactusProject()
                project.readXML(projPath)
            else:
                project = None
            row = params.asRow()
            row.extend(self.__jtStats(jtXmlRoot))
            row.extend(self.__totalAggregate(mafXmlRoot))
            row.extend(self.__growthStats(project))
            row.extend(self.__cactusTreeStats(tsXmlRoot))
            row.extend(self.__speciesAggregate(mafXmlRoot))
            rowstring = str(row)
            self.table.append(row)
    
    def addEmptyLine(self):
        row = self.getRows().next()
        emptyRow = []
        for i in xrange(len(row)):
            emptyRow.append("")
        self.table.append(emptyRow)
        
    def write(self, path):
        if len(self.table) > 0:
            emptyCols = self.__findEmptyColumns()
            outFile = open(path, "w")
            header = self.getHeader()
            outFile.write(self.__printLine(header, emptyCols))
            for row in self.getRows():
                assert len(row) == len(header)
                outFile.write(self.__printLine(row, emptyCols))
            outFile.close()
     
    def __printLine(self, rowAsList, emptyCols):
        rowString = ""
        for col in xrange(len(rowAsList)):
            if emptyCols[col] is False:
                rowString += "%s," % str(rowAsList[col])
        return "%s\n" % rowString[:len(rowString)-1] 
    
    def __pairTests(self, xmlRoot, i):
        assert i == 0 or i == 1
        homologyTest = xmlRoot.findall("homologyTests")[i]
        pairTests = homologyTest.find("homologyPairTests")
        for test in pairTests.findall("homologyTest"):
            yield test
    
    # returns [map of species names to aggregate sensitivity, 
    #          map of species names to aggregate specificity]
    def __speciesAggregate(self, xmlRoot):
        results = [dict(), dict()]
        for i in [0,1]:
            for test in self.__pairTests(xmlRoot, i):
                if test.attrib["sequenceB"] == "aggregate":
                    name = test.attrib["sequenceA"]
                    agg = test.find("aggregateResults").find("all")
                    val = agg.attrib["average"]
                    results[i][name] = val
                elif test.attrib["sequenceA"] == "aggregate":
                    name = test.attrib["sequenceB"]
                    agg = test.find("aggregateResults").find("all")
                    val = agg.attrib["average"]
                    results[i][name] = val
        return results            
        
    # returns [sensitivity, specificity]
    def __totalAggregate(self, xmlRoot):
        homologyTests = xmlRoot.findall("homologyTests")
        assert len(homologyTests) == 2
        results = []
        for ht in homologyTests:
            agg = ht.find("aggregateResults").find("all")
            results.append(float(agg.attrib["average"]))
        acc = (results[-1] + results[-2]) / 2
        results.append(acc)
        return results
    
    def __getElemAtt(self, node, elemName, attributeName, default=""):
            elem = node.find(elemName)
            if elem is not None:
                if attributeName in elem.attrib:
                    return elem.attrib[attributeName]
            return default
        
    # returns the total run time  
    def __jtStats(self, xmlRoot):       
        jtRow = [float(xmlRoot.attrib["total_run_time"]), 
                float(xmlRoot.attrib["total_clock"])]
        ttElem = xmlRoot.find("target_types")
        jtRow.append(self.__getElemAtt(ttElem, "CactusBaseLevelAlignerWrapper", "total_time"))
        jtRow.append(self.__getElemAtt(ttElem, "CactusBaseLevelAlignerWrapper", "total_clock"))
        jtRow.append(self.__getElemAtt(ttElem, "CactusBaseLevelAlignerWrapper", "max_clock"))
        jtRow.append(self.__getElemAtt(ttElem, "CactusCoreWrapper1", "total_time"))
        jtRow.append(self.__getElemAtt(ttElem, "CactusCoreWrapper1", "total_clock"))
        jtRow.append(self.__getElemAtt(ttElem, "CactusCoreWrapper1", "max_clock"))
        jtRow.append(self.__getElemAtt(ttElem, "CactusCoreWrapper2", "total_time"))
        jtRow.append(self.__getElemAtt(ttElem, "CactusCoreWrapper2", "total_clock"))
        jtRow.append(self.__getElemAtt(ttElem, "CactusCoreWrapper2", "max_clock"))
        jtRow.append(self.__getElemAtt(ttElem, "RunSelfBlast", "total_time"))
        jtRow.append(self.__getElemAtt(ttElem, "RunSelfBlast", "total_clock"))
        jtRow.append(self.__getElemAtt(ttElem, "RunSelfBlast", "max_clock"))
        jtRow.append(self.__getElemAtt(ttElem, "RunBlast", "total_time"))
        jtRow.append(self.__getElemAtt(ttElem, "RunBlast", "total_clock"))
        jtRow.append(self.__getElemAtt(ttElem, "RunBlast", "max_clock"))
        return jtRow
    
    # give some stats on how the reference geneomes grow:
    # root growth: ratio of root to biggest leaf
    # avg growth: average ratio of a genome's size to its children's
    def __growthStats(self, project):
        if project is None:
            return ["", ""]
        results = []
        tree = project.mcTree
        rootName = "Anc0" #tree.getRootName()
        rootExpPath = project.expMap[rootName]
        rootExp = ExperimentWrapper(ET.parse(rootExpPath).getroot())
        rootPath = rootExp.getReferencePath()
        rootSize = float(os.path.getsize(rootPath))
        leafNames = [tree.getName(i) for i in tree.getLeaves()]
        leafSizes = []
        for leaf in leafNames:
            leafPath = project.sequencePath(leaf)
            leafSize = float(os.path.getsize(leafPath))
            leafSizes.append(leafSize)
        results.append(rootSize / max(leafSizes))
        
        ratioSum = 0.0
        ratioCount = 0
        for expName, expPath in project.expMap.items():
            exp = ExperimentWrapper(ET.parse(expPath).getroot())
            rootPath = exp.getReferencePath()
            if os.path.exists(rootPath):
                rootSize = float(os.path.getsize(rootPath))
                for leafName, leafPath in exp.seqMap.items():
                    leafSize = float(os.path.getsize(leafPath))
                    ratio = rootSize / leafSize
                    ratioSum += ratio
                    ratioCount += 1
        avgRatio = ratioSum / ratioCount
        results.append(avgRatio)
        return results
    
    # get some stats on the cactus graph structure 
    # bbl = <chains><base_block_lengths >
    # tgs = <terminal_group_sizes>
    def __cactusTreeStats(self, treeStatsXmlRoot):
        if treeStatsXmlRoot is not None:
            results = []
            chains = treeStatsXmlRoot.find("chains")
            bbl = chains.find("base_block_lengths")
            results.append(bbl.attrib["min"])
            results.append(bbl.attrib["max"])
            results.append(bbl.attrib["avg"])
            
            tgs = treeStatsXmlRoot.find("terminal_group_sizes")
            results.append(tgs.attrib["min"])
            results.append(tgs.attrib["max"])
            results.append(tgs.attrib["avg"])
        else:
            results = ["NA", "NA", "NA", "NA", "NA", "NA"] 
        return results
                        
    # return boolean vector identifying empty columns
    def __findEmptyColumns(self):
        ec = []
        row = self.getRows().next()
        for i in row:
            ec.append(True)
        for row in self.getRows():
            colNum = 0
            for col in row:
                if col != "":
                    ec[colNum] = False
                colNum +=1
        return ec
    
    def getHeader(self):
        assert len(self.table) > 0
        row = self.table[0]
        header = copy.deepcopy(self.Header)
        species = []
        assert type(row[self.SensIdx] == 'dict')
        for entry, value in row[self.SensIdx].items():
            species.append(entry)
        species = sorted(species)
        for i in species:
            header.append("%s_sens" % i)
            header.append("%s_spec" % i)
        return header
    
    def getRows(self):
        for entry in self.table:
            row = entry[:self.SensIdx]
            species = []
            assert type(entry[self.SensIdx] == 'dict')
            assert type(entry[self.SpecIdx] == 'dict')
            for name, value in entry[self.SensIdx].items():
                species.append(name)
            species = sorted(species)
            for i in species:
                row.append(entry[self.SensIdx][i])
                row.append(entry[self.SpecIdx][i])
            yield row
    
def main():
    usage = "usage: %prog <mafcomp xml> <jobtree xml>"
    description = "TEST: print summary row for input"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--project", dest="projectPath", 
                      default = None, help="path to multi cactus project xml file")
    
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    summary = Summary()
    summary.addRow("name", Params(), args[1], args[0], options.projectPath)
    print summary.getHeader()
    for row in summary.getRows():
        print row
        

if __name__ == '__main__':    
    main()
    
    
