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

def getRootPathString():
    """
    function for finding external location
    """
    import progressiveBenchmarks.src.pipeline
    i = os.path.abspath(progressiveBenchmarks.src.pipeline.__file__)
    return os.path.split(os.path.split(i)[0])[0] #os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

# default parameters (ie read directly from config xml)
# vanilla is set to true once for each set of iteration parameters
class ParamsGenerator:
    def __init__(self):
        self.minChainLength = [None]
        self.minBlockDegree = [None]
        self.maxGroupSize = [None]
        self.outgroupStrategy = [None]
        self.singleCopyStrategy = [None]
        self.requiredFraction = [None]
        self.selfAlignment = [None]
        self.subtreeSize = [None]
        self.vanilla = [True, False]
        self.kyotoTycoon = [None]
        self.templatePath = [None]
        
    def generate(self):
        for tp in self.templatePath:   
            for mc in self.minChainLength:
                for mb in self.minBlockDegree:
                    for mg in self.maxGroupSize:
                        for va in self.vanilla:
                            if va == True:
                                params = Params()
                                params.templatePath = tp
                                params.minChainLength = mc
                                params.minBlockDegree = mb
                                params.maxGroupSize = mg
                                params.vanilla = va
                                yield params
                            else:     
                                for og in self.outgroupStrategy:
                                    for sc in self.singleCopyStrategy:
                                        for cf in self.requiredFraction:
                                            for sa in self.selfAlignment:
                                                for st in self.subtreeSize:
                                                    for kt in self.kyotoTycoon:
                                                        params = Params()
                                                        params.templatePath = tp
                                                        params.minChainLength = mc
                                                        params.minBlockDegree = mb
                                                        params.maxGroupSize = mg
                                                        params.outgroupStrategy = og
                                                        params.singleCopyStrategy = sc
                                                        params.subtreeSize = st
                                                        params.requiredFraction = cf
                                                        params.selfAlignment = sa
                                                        params.vanilla = va
                                                        params.kyotoTycoon = kt
                                                        yield params
                                
class EverythingButSelf():
    class EverythingButSelf_MB2(ParamsGenerator):
        def __init__(self):
            ParamsGenerator.__init__(self)
            self.minChainLength = [32, 64, 128, 256, 512, 1024]
            self.minBlockDegree = [2]
            self.maxGroupSize = [10000000000]
            self.outgroupStrategy = ['none', 'greedy', 'greedyLeaves']
            self.singleCopyStrategy = ['none', 'outgroup', 'all']
            self.requiredFraction = [0, 0.67, 1]

    class EverythingButSelf_MB0(ParamsGenerator):
        def __init__(self):
            ParamsGenerator.__init__(self)
            self.minChainLength = [32, 64, 128, 256, 512, 1024]
            self.minBlockDegree = [0]
            self.maxGroupSize = [10000000000]
            self.outgroupStrategy = ['none']
            self.singleCopyStrategy = ['none', 'all']
            self.requiredFraction = [0]
   
    def generate(self):
        for p in self.EverythingButSelf_MB2().generate():
            yield p
        for p in self.EverythingButSelf_MB0().generate():
            yield p
        
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
        
class SingleCase(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.minChainLength = [ 256 ]
        self.outgroupStrategy = ['greedy']
        self.singleCopyStrategy = ['outgroup']
        self.requiredFraction = [0.67]
        self.selfAlignment = [False]
        self.vanilla = [ False ]
        self.minBlockDegree = [ 2]

class KyotoTycoon(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.kyotoTycoon = [True]
        self.outgroupStrategy = ['greedy', 'greedyLeaves']
        self.singleCopyStrategy = ['outgroup']
        self.selfAlignment = [False]
        self.vanilla = [False]
        
class LastzTuning(ParamsGenerator):
    def lastzPath(self, x):
        return os.path.join(getRootPathString(), "lib", "lastz_tuning", 
                            "config_lastz_%d.xml" % x) 
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.templatePath = [self.lastzPath(1)]
        self.vanilla = [True, False]
                             
    