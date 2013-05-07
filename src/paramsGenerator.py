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
        self.annealingRounds = [None]
        self.minBlockDegree = [None]
        self.repeatMask = [None]
        self.outgroupStrategy = [None]
        self.outgroupThreshold = [None]
        self.selfAlignment = [None]
        self.subtreeSize = [None]
        self.vanilla = [True, False]
        self.kyotoTycoon = [None]
        self.templatePath = [None]
        self.numThreads = [None]
        
    def generate(self):
        for tp in self.templatePath:   
            for ar in self.annealingRounds:
                for mb in self.minBlockDegree:
                    for rm in self.repeatMask:
                        for nt in self.numThreads:
                            for va in self.vanilla:
                                if va == True:
                                    params = Params()
                                    params.templatePath = tp
                                    params.annealingRounds = ar
                                    params.minBlockDegree = mb
                                    params.repeatMask = rm
                                    params.numThreads = nt
                                    params.vanilla = va
                                    yield params
                                else:     
                                    for og in self.outgroupStrategy:
                                        for ot in self.outgroupThreshold:
                                            if og != 'greedy' and self.outgroupThreshold.index(ot) > 0:
                                                continue
                                            for sa in self.selfAlignment:
                                                for st in self.subtreeSize:
                                                    for kt in self.kyotoTycoon:
                                                        params = Params()
                                                        params.templatePath = tp
                                                        params.annealingRounds = ar
                                                        params.minBlockDegree = mb
                                                        params.repeatMask = rm
                                                        params.outgroupStrategy = og
                                                        params.outgroupThreshold = ot
                                                        params.subtreeSize = st
                                                        params.selfAlignment = sa
                                                        params.vanilla = va
                                                        params.kyotoTycoon = kt
                                                        params.numThreads = nt
                                                        yield params
                        
class EverythingButSelf():
    class EverythingButSelf_MB2(ParamsGenerator):
        def __init__(self):
            ParamsGenerator.__init__(self)
            self.minBlockDegree = [2]
            self.outgroupStrategy = ['none', 'greedy', 'greedyLeaves']

    class EverythingButSelf_MB0(ParamsGenerator):
        def __init__(self):
            ParamsGenerator.__init__(self)
            self.minBlockDegree = [0]
            self.outgroupStrategy = ['none']
   
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
        self.selfAlignment = [True, False]

class BasicProgressive(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.outgroupStrategy = ['none', 'greedy']
        self.selfAlignment = [True, False]
    
class SmallProgressive(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.outgroupStrategy = ['none', 'greedy']
        
class SingleCase(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.outgroupStrategy = ['greedy']
        #self.outgroupThreshold = [0]
        self.selfAlignment = [False]
        self.vanilla = [ False ]
        self.repeatMask = [10]
        #self.minBlockDegree = [ 2 ]

class KyotoTycoon(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.kyotoTycoon = [ "inMemory" ] #[True, "inMemory", "inMemoryNoSnapshot"]
        self.outgroupStrategy = ['greedy' ]
        self.outgroupThreshold = [0]
        self.selfAlignment = [False]
        self.vanilla = [False]
        self.repeatMask = [10]
        
class RepeatMasking(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.kyotoTycoon = [ "inMemoryNoSnapshot" ] #[True, "inMemory", "inMemoryNoSnapshot"]
        self.outgroupStrategy = ['greedyLeaves' ]
        self.outgroupThreshold = [0]
        self.selfAlignment = [False]
        self.vanilla = [False]
        self.repeatMask = [ 3, 5, 10 ]
        
class LastzTuning(ParamsGenerator):
    def lastzPath(self, x):
        return os.path.join(getRootPathString(), "lib", "lastz_tuning", 
                            "config_lastz_%d.xml" % x) 
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.repeatMask = [None, 10]
        #self.repeatMask = [10]
        self.templatePath = [self.lastzPath(i) for i in range(1,6)]
        #self.templatePath = [self.lastzPath(1)]
        self.annealingRounds = ["2 3 4 8 16 32 64 128", "2 3 4 8 16 32 64 128 256 512", "8 128 512"]
        #self.annealingRounds = ["2 3 4 8 16 32 64 128"]
        self.vanilla = [False]

class DefaultProgTest(ParamsGenerator):
    def __init__(self):
        ParamsGenerator.__init__(self)
        self.vanilla = [False]
        
class NumThreadsTest(SingleCase):
    def __init__(self):
        SingleCase.__init__(self)
        self.numThreads = [1, 2, 3, 4]
        self.outgroupStrategy = ['greedy', 'none']
                             
    