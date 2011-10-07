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
    
    def __init__(self):
        self.minChainLength = None
        self.minBlockDegree = None
        self.maxGroupSize = None
        self.iterationParams = None
        self.outgroupStrategy = None
        self.singleCopyStrategy = None
        self.requiredFraction = None
        self.selfAlignment = None
        self.subtreeSize = None
        self.vanilla = False
    
    # only write non-None attributes, idea being
    # that defaults are already in the config.        
    def applyToXml(self, config):
        def setAtt(elem, atName, val):
            if val is not None:
                elem.attrib[atName] = str(val)
        
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
        iterationElem = iterationsElem.findall("iteration")[-1]        
        setAtt(iterationElem, "minimumChainLength", self.minChainLength)
        setAtt(iterationElem, "minimumBlockDegree", self.minBlockDegree)
        setAtt(iterationElem, "maximumGroupSize", self.maxGroupSize)                    

    def check(self):
        if self.doVanilla == True:
            assert self.outgroupStrategy is None
            assert self.singleCopyStrategy is None
            assert self.requiredFraction is None
            assert self.selfAlignment is None
            assert self.subtreeSize is None
            
    def __str__( self ):
        def printItem(name, value):
            if value is None:
                return ''
            else:
                return "_%s%s" % (name, str(value).title())
        
        self.check()    
        token = ""
        
        token += printItem("mc", self.minChainLength)
        token += printItem("mb", self.minBlockDegree)
        token += printItem("mg", self.maxGroupSize)
        token += printItem("og", self.outgroupStrategy)
        token += printItem("sc", self.singleCopyStrategy)
        token += printItem("cf", self.requiredFraction)
        token += printItem("sa", self.selfAlignment)
        token += printItem("st", self.subtreeSize)
        if self.doVanilla:
            token += printItem("", "vanilla")
        if token == "":
            token = "_Default"
        return token
            
                