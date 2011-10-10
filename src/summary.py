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

class Summary:
    Header = Params.Header + \
    ["Run_Time", "Clock_Time", "Sensitivity", "Specificity"]
    SensIdx = len(Header)
    SpecIdx = SensIdx + 1 
    def __init__(self):
        self.table = []

    def addRow(self, catName, params, jobTreeStatsPath, mafCompPath):
        if os.path.isfile(jobTreeStatsPath) and os.path.isfile(mafCompPath):
            mafXmlRoot = ET.parse(mafCompPath).getroot()
            jtXmlRoot = ET.parse(jobTreeStatsPath).getroot()
            row = params.asRow()
            row.extend(self.__jtStats(jtXmlRoot))
            row.extend(self.__totalAggregate(mafXmlRoot))
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
        homologyTest = xmlRoot.findall("homology_tests")[i]
        pairTests = homologyTest.find("homology_pair_tests")
        for test in pairTests.findall("homology_test"):
            yield test
    
    # returns [map of species names to aggregate sensitivity, 
    #          map of species names to aggregate specificity]
    def __speciesAggregate(self, xmlRoot):
        results = [dict(), dict()]
        for i in [0,1]:
            for test in self.__pairTests(xmlRoot, i):
                if test.attrib["sequenceB"] == "aggregate":
                    name = test.attrib["sequenceA"]
                    agg = test.find("aggregate_results").find("all")
                    val = agg.attrib["average"]
                    results[i][name] = val
                elif test.attrib["sequenceA"] == "aggregate":
                    name = test.attrib["sequenceB"]
                    agg = test.find("aggregate_results").find("all")
                    val = agg.attrib["average"]
                    results[i][name] = val
        return results            
        
    # returns [sensitivity, specificity]
    def __totalAggregate(self, xmlRoot):
        homologyTests = xmlRoot.findall("homology_tests")
        assert len(homologyTests) == 2
        results = []
        for ht in homologyTests:
            agg = ht.find("aggregate_results").find("all")
            results.append(float(agg.attrib["average"]))
        return results
    
    # returns the total run time  
    def __jtStats(self, xmlRoot):
        return [float(xmlRoot.attrib["total_run_time"]), 
                float(xmlRoot.attrib["total_clock"])]
    
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
    
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    summary = Summary()
    summary.addRow("name", "options", args[1], args[0])
    print summary.getHeader()
    for row in summary.getRows():
        print row
        

if __name__ == '__main__':    
    main()
    
    