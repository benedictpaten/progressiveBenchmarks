#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" generate sets of test parameters

"""

import os
import xml.etree.ElementTree as ET
import sys

from progressiveBenchmarks.src.params import Params

# default parameters (ie read directly from config xml)
# doVanilla is set to true once for each set of iteration parameters
class ParamsGenerator:
    def __init__(self):
        self.iterationParams = [None]
        self.outgroupStrategy = [None]
        self.singleCopyStrategy = [None]
        self.requiredFraction = [None]
        self.selfAlignment = [None]
        self.subtreeSize = [None]
        self.doVanilla = [True, False]
        
    def generate(self):   
        for it in self.iterationParams:
            for va in self.doVanilla:
                if va == True:
                    params = Params()
                    params.iterationParams = it
                    params.doVanilla = va
                    yield params
                else:     
                    for og in self.outgroupStrategy:
                        for sc in self.singleCopyStrategy:
                            for cf in self.requiredFraction:
                                for sa in self.selfAlignment:
                                    for st in self.subtreeSize:
                                        params = Params()
                                        params.iterationParams = it
                                        params.outgroupStrategy = og
                                        params.singleCopyStrategy = sc
                                        params.subtreeSize = st
                                        params.requiredFraction = cf
                                        params.selfAlignment = sa
                                        params.doVanilla = va
                                        yield params
                

# all the progressive-related combinations                                
class AllProgressive(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.outgroupStrategy = ['none', 'greedy', 'greedyLeaves']
        self.singleCopyStrategy = ['none', 'outgroup', 'all']
        self.requiredFraction = [0, 0.67, 1]
        self.selfAlignment = [True, False]

class BasicProgressive(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.outgroupStrategy = ['none', 'greedy']
        self.singleCopyStrategy = ['outgroup']
        self.requiredFraction = [0]
        self.selfAlignment = [True, False]
    
class SmallProgressive(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.outgroupStrategy = ['none', 'greedy']

    