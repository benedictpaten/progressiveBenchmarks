#!/usr/bin/env python

"""Script runs cactus to build a bunch of reference genome assemblies for a given locus
"""

import os
import xml.etree.ElementTree as ET
import xml
import sys
from optparse import OptionParser

from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import setLoggingFromOptions

from cactus.shared.config import CactusWorkflowExperiment
from cactus.shared.common import runCactusWorkflow

from cactusTools.shared.common import runCactusTreeStats
from cactusTools.shared.common import runCactusMAFGenerator

from sonLib.bioio import getTempFile, getTempDirectory
from sonLib.bioio import fastaRead, fastaWrite
from sonLib.bioio import system

from cactus.shared.test import getCactusInputs_blanchette

from jobTree.test.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

def getRootPathString():
    """
    function for finding external location
    """
    import os
    import referenceScripts.bin.pipeline
    i = os.path.abspath(referenceScripts.bin.pipeline.__file__)
    return os.path.split(os.path.split(i)[0])[0] #os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

def getCactusDiskString(alignmentFile):
    return "<st_kv_database_conf type=\"tokyo_cabinet\"><tokyo_cabinet database_dir=\"%s\"/></st_kv_database_conf>" % alignmentFile

class MakeAlignment(Target):
    """Target runs the alignment.
    """
    def __init__(self, options,
                 sequences, 
                 newickTree,
                 outputDir, 
                 useOutgroup,
                 doSelfAlignment, 
                 trueAlignment=None):
                 #requiredSpecies,
                 #singleCopySpecies,
                 #referenceAlgorithm, minimumBlockDegree, 
                 #blastAlignmentString, baseLevel, maxNumberOfChains, permutations,
                 #theta, useSimulatedAnnealing, heldOutSequence):
        Target.__init__(self, cpu=1, memory=4000000000)
        self.sequences = sequences
        self.newickTree = newicktree,
        #self.requiredSpecies = requiredSpecies
        #self.singleCopySpecies = singleCopySpecies
        self.outputDir = outputDir
        #self.referenceAlgorithm = referenceAlgorithm
        #self.minimumBlockDegree = int(minimumBlockDegree)
        #self.blastAlignmentString = blastAlignmentString
        #self.baseLevel = baseLevel
        #self.maxNumberOfChains = maxNumberOfChains
        #self.permutations = permutations
        #self.theta = theta
        #self.useSimulatedAnnealing = useSimulatedAnnealing
        self.options = options
        #self.heldOutSequence = heldOutSequence
        self.useOutgroup = useOutgroup
        self.doSelfAlignment = doSelfAlignment
        self.trueAlignment = trueAlignment
    
    def run(self):
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)
        cactusAlignmentName = "cactusAlignment"
        outputFile = os.path.join(self.outputDir, cactusAlignmentName)
        if not os.path.exists(outputFile):
            config = ET.parse(os.path.join(getRootPathString(), "lib", "cactus_workflow_config.xml")).getroot()
            
            #Set the parameters
            
            
            #Write the config file
            tempConfigFile = os.path.join(self.getLocalTempDir(), "config.xml")
            fileHandle = open(tempConfigFile, 'w')
            tree = ET.ElementTree(config)
            tree.write(fileHandle)
            fileHandle.close()
         
            #Make the experiment file
            tempExperimentFile = os.path.join(self.getLocalTempDir(), "experiment.xml")
            cactusWorkflowExperiment = CactusWorkflowExperiment(
                                                 sequences=self.sequences.split(), 
                                                 newickTreeString=self.newickTree, 
                                                 #requiredSpecies=self.requiredSpecies,
                                                 #singleCopySpecies=self.singleCopySpecies,
                                                 databaseName=cactusAlignmentName,
                                                 outputDir=self.getLocalTempDir(),
                                                 configFile=tempConfigFile)
            cactusWorkflowExperiment.writeExperimentFile(tempExperimentFile)
            
            #The jobtree
            tempJobTreeDir = os.path.join(self.getLocalTempDir(), "jobTree")
            
            #The place to put the temporary experiment dir
            tempDir = getTempDirectory(os.getcwd())
            tempExperimentDir = os.path.join(tempDir, "exp")
            
            #The temporary experiment 
            runCactusCreateMultiCactusProject(tempExperimentFile, 
                                              tempExperimentDir,
                                              useOutgroup=self.useOutgroup,
                                              doSelfAlignment=self.doSelfAlignment)
            logger.info("Setup the cactus progressive experiment")
            
            runCactusProgressive(os.path.join(tempExperimentDir, "exp_project.xml"), 
                                 tempJobTreeDir, 
                                 batchSystem=batchSystem, 
                                 buildMAF=True,
                                 joinMAF=True,
                                 #buildTrees=buildTrees, buildFaces=buildFaces, buildReference=buildReference,
                                 jobTreeStats=True)
            logger.info("Ran the progressive workflow")
            
            #Check if the jobtree completed sucessively.
            runJobTreeStatusAndFailIfNotComplete(tempJobTreeDir)
            logger.info("Checked the job tree dir")
            
            #Now copy the true assembly back to the output
            system("mv %s %s/experiment.xml" % (tempExperimentFile, self.outputDir))
            system("mv %s %s/config.xml" % (tempConfigFile, self.outputDir))
            system("mv %s %s" % (tempExperimentDir, self.outputDir))
            system("jobTreeStats --jobTree %s --outputFile %s/jobTreeStats.xml" % (tempJobTreeDir, self.outputDir))
            #We're done!
        #self.addChildTarget(MakeStats(outputFile, self.outputDir, self.options)) 

class MakeBlanchetteAlignments(Target):
    def __init__(self, options, useOutgroup, doSelfAlignment):
        self.options = options
        self.useOutgroup = useOutgroup
        self.doSelfAlignment = doSelfAlignment
    
    def run(self):
        outputDir = "blanchette-%s-%s" % (self.useOutgroup, self.doSelfAlignment)
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)
        for i in xrange(5):
            sequences, newickTreeString = getCactusInputs_blanchette(i)
            outputDir = os.path.join(outputDir, str(i))
            self.addChild(MakeAlignment(self.options, sequences, newickTree, outputDir, 
                                        self.useOutgroup, self.doSelfAlignment))
        #self.setFollowOnJob(MakeStats(options, outputDir))
        
class MakeBlanchetteStats(MakeAlignment):
    """Builds basic stats and the maf alignment.
    """
    def run(self):
        previousOutputFile = None
        blanchettePath = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation")
        for i in xrange(len(self.sequences)):
            trueAlignmentMFA = os.path.join(os.path.join(blanchettePath, "%.2i.job" % i), "true.mfa")
            trueAlignmentMAF = os.path.join(self.getLocalTempDir(), "temp.maf")
            system("mfaToMaf --mfaFile %s --outputFile %s" % (trueAlignmentMFA, trueAlignmentMAF))
            predictedAlignment = os.path.join(os.path.join(self.outputDir, str(i)), "alignment.maf")
            outputFile = os.path.join(self.getLocalTempDir(), "temp%i" % i)
            system("mafComparator --mafFile1 %s --mafFile2 %s --outputFile %s" % (trueAlignment, predictedAlignment, outputFile))
            if previousOutputFile != None:
                system("mergeMafComparatorResults.py --results1 %s --results2 %s --outputFile %s" % (outputFile, previousOutputFile, outputFile))
            previousOutputFile = outputFile
        system("mv %s %s" % (previousOutputFile, os.path.join(self.outputDir, "mafComparison.xml")))   
      
class MakeAllAlignments(Target):
    """Makes alignments using pipeline.
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options
    
    def run(self):
        for useOutgroup in (True, False):
            for doSelfAlignment in (True, False):
                self.addChildTarge(MakeBlanchetteAlignments(self.options, useOutgroup, doSelfAlignment))
                #self.addChildTarget(MakeEvolverPrimates(self.options, useOutgroup, doSelfAlignment))

def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
    parser.add_option("--outputDir", dest="outputDir")
    
    Stack.addJobTreeOptions(parser)
    
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))
    
    Stack(MakeAllAlignments(options)).startJobTree(options)
    logger.info("Done with job tree")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from referenceScripts.bin.pipeline import *
    _test()
    main()

