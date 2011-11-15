#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" store the different configuration parameters for
progressive cactus

"""

import os
import xml.etree.ElementTree as ET
import sys
import copy

class Params:
    
    Header = ["Style", "Template", "MinChainLen", "MinBlockDeg", "MaxGroupSize", \
              "Outgroup", "SingleCpy", "ReqFrac", "Self", "SubtreeSize", "Kyoto"]
                
    def __init__(self):
        self.minChainLength = None
        self.minBlockDegree = None
        self.maxGroupSize = None
        self.outgroupStrategy = None
        self.singleCopyStrategy = None
        self.requiredFraction = None
        self.selfAlignment = None
        self.subtreeSize = None
        self.vanilla = False
        self.kyotoTycoon = None
        self.templatePath = None
    
    # only write non-None attributes, idea being
    # that defaults are already in the config.        
    def applyToXml(self, xmlTree):
        def setAtt(elem, atName, val):
            if val is not None:
                elem.attrib[atName] = str(val)
        
        if self.templatePath != None:
            xmlTree._setroot(ET.parse(self.templatePath).getroot())
           
        config = xmlTree.getroot()
                
        mcElem = config.find("multi_cactus")
        ogElem = mcElem.find("outgroup")
        setAtt(ogElem, "strategy", self.outgroupStrategy)
        coverageElem = mcElem.find("coverage")
        setAtt(coverageElem, "required_fraction", self.requiredFraction)
        setAtt(coverageElem, "single_copy_strategy", self.singleCopyStrategy)
        decompElem = mcElem.find("decomposition")
        setAtt(decompElem, "self_alignment", self.selfAlignment)
        setAtt(decompElem, "subtree_size", self.subtreeSize)
    
        iterationsElem = config.find("alignment").find("iterations")
        iterations = iterationsElem.findall("iteration")
        cafElem = iterations[-2]
        assert cafElem.attrib["type"] == "blast"
        assert cafElem.attrib["number"] == "0"
        coreElem = cafElem.find("core")
        setAtt(coreElem, "minimumChainLength", self.minChainLength)
        setAtt(coreElem, "maximumGroupSize", self.maxGroupSize)
        
        barElem = iterations[-1]
        assert barElem.attrib["type"] == "base"
        assert barElem.attrib["number"] == "1"
        setAtt(barElem, "minimumBlockDegree", self.minBlockDegree)                    

    def check(self):
        if self.vanilla == True:
            assert self.outgroupStrategy is None
            assert self.singleCopyStrategy is None
            assert self.requiredFraction is None
            assert self.selfAlignment is None
            assert self.subtreeSize is None
            assert self.kyotoTycoon is None or self.kyotoTycoon is False
            
    def __str__( self ):
        def printItem(name, value):
            if value is None:
                return ''
            else:
                return "_%s%s" % (name, str(value).title())
        
        self.check()    
        token = ""
        
        tpName = self.templatePath
        if tpName is not None:
            tpName = os.path.basename(tpName)
            tpName = os.path.splitext(tpName)[0]
        
        token += printItem("tp", tpName)
        token += printItem("mc", self.minChainLength)
        token += printItem("mb", self.minBlockDegree)
        token += printItem("mg", self.maxGroupSize)
        token += printItem("og", self.outgroupStrategy)
        token += printItem("sc", self.singleCopyStrategy)
        token += printItem("cf", self.requiredFraction)
        token += printItem("sa", self.selfAlignment)
        token += printItem("st", self.subtreeSize)
        token += printItem("kt", self.kyotoTycoon)
        if self.vanilla:
            token += printItem("", "vanilla")
        if token == "":
            token = "_Default"
        return token
    
    def asRow(self):
        def addItem(row, value):
            if value is None:
                row.append("")
            else:
                row.append(value)
            
        row = []
        name = "Progressive"
        if self.vanilla is True:
            name = "Vanilla"
            
        tpName = self.templatePath
        if tpName is not None:
            tpName = os.path.basename(tpName)
            tpName = os.path.splitext(tpName)[0]
            
        addItem(row, name) 
        addItem(row, tpName)   
        addItem(row, self.minChainLength)
        addItem(row, self.minBlockDegree)
        addItem(row, self.maxGroupSize)
        addItem(row, self.outgroupStrategy)
        addItem(row, self.singleCopyStrategy)
        addItem(row, self.requiredFraction)
        addItem(row, self.selfAlignment)
        addItem(row, self.subtreeSize)
        addItem(row, self.kyotoTycoon)
        return row

        
                