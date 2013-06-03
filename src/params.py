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
    
    Header = ["Style", "Template", "AnnealingRounds", "MinBlockDeg", 
              "RepeatMask", "Outgroup", "OgThreshold", "Self", \
              "SubtreeSize", "Kyoto", "numThreads"]
                
    def __init__(self):
        self.annealingRounds = None
        self.minBlockDegree = None
        self.repeatMask = None
        self.outgroupStrategy = None
        self.outgroupThreshold = None
        self.selfAlignment = None
        self.subtreeSize = None
        self.vanilla = False
        self.kyotoTycoon = None
        self.templatePath = None
        self.numThreads = None
    
    def setAtt(self, elem, atName, val):
        if val is not None:
            elem.attrib[atName] = str(val)

    # only write non-None attributes, idea being
    # that defaults are already in the config.        
    def applyToXml(self, xmlTree):        
        if self.templatePath != None:
            xmlTree._setroot(ET.parse(self.templatePath).getroot())
           
        config = xmlTree.getroot()
                
        mcElem = config.find("multi_cactus")
        ogElem = mcElem.find("outgroup")
        self.setAtt(ogElem, "strategy", self.outgroupStrategy)
        self.setAtt(ogElem, "threshold", self.outgroupThreshold)
        decompElem = mcElem.find("decomposition")
        self.setAtt(decompElem, "self_alignment", self.selfAlignment)
        self.setAtt(decompElem, "subtree_size", self.subtreeSize)
    
        cafElem = config.find("caf")
        if self.annealingRounds is not None:
            self.setAtt(cafElem, "annealingRounds", self.annealingRounds)
            trim = cafElem.attrib["trim"].split()[0]
            trimList = [trim] * len(self.annealingRounds.split())
            self.setAtt(cafElem, "trim", " ".join(trimList))
               
        barElem = config.find("bar") 
        self.setAtt(barElem, "minimumBlockDegree", self.minBlockDegree)
        self.setAtt(barElem, "numThreads", self.numThreads)
        
    def check(self):
        if self.vanilla == True:
            assert self.outgroupStrategy is None
            assert self.outgroupThreshold is None
            assert self.selfAlignment is None
            assert self.subtreeSize is None
            assert self.kyotoTycoon is None or self.kyotoTycoon is False
            
    def __str__( self ):
        def printItem(name, value):
            if value is None:
                return ''
            else:
                return ("_%s%s" % (name, str(value).title())).replace(" ", ".")
        
        self.check()    
        token = ""
        
        tpName = self.templatePath
        if tpName is not None:
            tpName = os.path.basename(tpName)
            tpName = os.path.splitext(tpName)[0]
        
        token += printItem("tp", tpName)
        token += printItem("ar", self.annealingRounds)
        token += printItem("mb", self.minBlockDegree)
        token += printItem("rm", self.repeatMask)
        token += printItem("og", self.outgroupStrategy)
        token += printItem("ot", self.outgroupThreshold)
        token += printItem("sa", self.selfAlignment)
        token += printItem("st", self.subtreeSize)
        token += printItem("kt", self.kyotoTycoon)
        token += printItem("nt", self.numThreads)
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
        addItem(row, self.annealingRounds)
        addItem(row, self.minBlockDegree)
        addItem(row, self.repeatMask)
        addItem(row, self.outgroupStrategy)
        addItem(row, self.outgroupThreshold)
        addItem(row, self.selfAlignment)
        addItem(row, self.subtreeSize)
        addItem(row, self.kyotoTycoon)
        addItem(row, self.numThreads)
        return row

        
                
