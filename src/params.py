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


class Params:
    class IterationParams:
        def __init__(self):
            self.minChainLength = None
            self.minBlockDegree = None
            self.maxGroupSize = None
            
    def __init__(self):
        self.iterationParams = None
        self.outgroupStrategy = None
        self.singleCopyStrategy = None
        self.requiredFraction = None
        self.selfAlignment = None
        self.subtreeSize = None
        self.doVanilla = False
    
    # only write non-None attributes, idea being
    # that defaults are already in the config.        
    def applyToXml(self, config):
        def setAtt(elem, atName, val):
            if val is not None:
                elem.attrib[atName] = str(val)
    
        mcElem = config.find("multi_cactus")
        ogElem = mcElem.find("outgroup")
        ogElem.attrib["strategy"] = self.outgroupStrategy
        coverageElem = mcElem.find("coverage")
        setAtt(coverageElem, "required_fraction", self.requiredFraction)
        setAtt(coverageElem, "single_copy_strategy", self.singleCopyStrategy)
        decompElem = mcElem.find("decomposition")
        setAtt(decompElem, "self_alignment", self.selfAlignment)
        setAtt(decompElem, "subtree_size", self.subtreeSize)
        
        if self.iterationParams is not None:
            iterationsElem = config.find("alignment").find("iterations")
            iterationNumber = 0
            for iterationElem in iterationsElem.findall("iteration"):
                if iterationNumber < len(self.iterationParams):
                    setAtt(iterationElem, "minimumChainLength", 
                           self.iterationParams[iterationNumber].minChainLength)
                    setAtt(iterationElem, "minimumBlockDegree",
                           self.iterationParams[iterationNumber].minBlockDegree)
                    setAtt(iterationElem, "maximumGroupSize",
                           self.iterationParams[iterationNumber].maxGroupSize)                   
                iterationNumber += 1
        
    def __str__( self ):
        def printItem(name, value):
            if value is None:
                return ''
            else:
                return "_%s%s" % (name, str(value))
        
        token = ""
        if self.iterationParams is not None:
            token += printItem("mc", self.iterationParams.minChainLength)
            token += printItem("mb", self.iterationParams.minBlockDegree)
            token += printItem("mg", self.iterationParams.maxGroupSize)
        token += printItem("og", self.outgroupStrategy)
        token += printItem("sc", self.singleCopyStrategy)
        token += printItem("cf", self.requiredFraction)
        token += printItem("sa", self.selfAlignment)
        token += printItem("st", self.subtreeSize)
        token += printItem("dv", self.doVanilla)
        return token
            
                