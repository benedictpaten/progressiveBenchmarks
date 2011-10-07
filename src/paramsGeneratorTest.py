#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
"""
"""

import unittest
import os
import sys
import copy
import xml.etree.ElementTree as ET
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system

from progressiveBenchmarks.src.paramsGenerator import AllProgressive

class TestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def testAllProgressiveUnique(self):
        names = set()
        for params in AllProgressive().generate():
            name = str(params)
            print name
            assert name not in names
    
    def testAllDoVanillaOnce(self):
        names = set()
        vinCount = 0
        for params in AllProgressive().generate():
            if params.doVanilla == True:
                vinCount += 1
        print vinCount
        assert vinCount == 1
            
def main():
    unittest.main()
    
if __name__ == '__main__':
    main()
