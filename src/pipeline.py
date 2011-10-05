#!/usr/bin/env python

"""Script runs cactus to build a bunch of reference genome assemblies for a given locus
"""

import os
import xml.etree.ElementTree as ET
import xml
import sys
import re
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

from cactus.shared.common import runCactusProgressive
from cactus.shared.common import runCactusCreateMultiCactusProject
from cactus.shared.test import getInputs

from sonLib.bioio import TestStatus

from progressiveBenchmarks.src.params import Params
from progressiveBenchmarks.src.paramsGenerator import AllProgressiveCombos
from progressiveBenchmarks.src.paramsGenerator import SmallProgressive
from progressiveBenchmarks.src.applyNamingToMaf import applyNamingToMaf

def getRootPathString():
    """
    function for finding external location
    """
    import progressiveBenchmarks.src.pipeline
    i = os.path.abspath(progressiveBenchmarks.src.pipeline.__file__)
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
                 params):
                 #requiredSpecies,
                 #singleCopySpecies,
                 #referenceAlgorithm, minimumBlockDegree, 
                 #blastAlignmentString, baseLevel, maxNumberOfChains, permutations,
                 #theta, useSimulatedAnnealing, heldOutSequence):
        Target.__init__(self, cpu=4, memory=4000000000)
        self.sequences = sequences
        self.newickTree = newickTree
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
        self.params = params
    
    def run(self):
        logger.debug("Going to put the alignment in %s" % self.outputDir)
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)

        if not os.path.exists(os.path.join(self.outputDir, "progressiveCactusAlignment")):
            config = ET.parse(os.path.join(getRootPathString(), "lib", "cactus_workflow_config.xml")).getroot()
            
            #Set the parameters
            tempLocalDir = os.path.join(self.outputDir, "tempProgressiveCactusAlignment")
            system("rm -rf %s" % tempLocalDir)
            os.mkdir(tempLocalDir)
            
            #Set the config parameters
            self.params.applyToXml(config)
        
            #Write the config file
            tempConfigFile = os.path.join(tempLocalDir, "config.xml")
            fileHandle = open(tempConfigFile, 'w')
            tree = ET.ElementTree(config)
            tree.write(fileHandle)
            fileHandle.close()
         
            #Make the experiment file
            tempExperimentFile = os.path.join(tempLocalDir, "experiment.xml")
            
            cactusWorkflowExperiment = CactusWorkflowExperiment(
                                                 sequences=self.sequences, 
                                                 newickTreeString=self.newickTree, 
                                                 #requiredSpecies=self.requiredSpecies,
                                                 #singleCopySpecies=self.singleCopySpecies,
                                                 databaseName="cactusAlignment",
                                                 outputDir=tempLocalDir,
                                                 configFile=tempConfigFile)
            cactusWorkflowExperiment.writeExperimentFile(tempExperimentFile)
            
            #The jobtree
            tempJobTreeDir = os.path.join(tempLocalDir, "jobTree")
            
            #The place to put the temporary experiment dir
            tempExperimentDir = os.path.join(tempLocalDir, "progressiveCactusAlignment")
            
            #The temporary experiment 
            runCactusCreateMultiCactusProject(tempExperimentFile, 
                                              tempExperimentDir)
            logger.info("Setup the cactus progressive experiment")
            
            runCactusProgressive(os.path.join(tempExperimentDir, "progressiveCactusAlignment_project.xml"), 
                                 tempJobTreeDir, 
                                 #batchSystem=batchSystem, 
                                 buildMaf=True,
                                 joinMaf=True,
                                 #buildTrees=buildTrees, buildFaces=buildFaces, buildReference=buildReference,
                                 jobTreeStats=True,
                                 maxThreads=4)
            logger.info("Ran the progressive workflow")
            
            #Check if the jobtree completed sucessively.
            runJobTreeStatusAndFailIfNotComplete(tempJobTreeDir)
            logger.info("Checked the job tree dir for the progressive run")
            
            #Now copy the true assembly back to the output
            system("mv %s %s/experiment.xml" % (tempExperimentFile, self.outputDir))
            system("mv %s %s" % (tempExperimentDir, self.outputDir))
            system("jobTreeStats --jobTree %s --outputFile %s/jobTreeStats.xml" % (tempJobTreeDir, self.outputDir))
                
            if self.params.doVanilla:
                #Now do standard cactus..
                #Make the experiment file
                tempExperimentFile2 = os.path.join(tempLocalDir, "experimentVanilla.xml")

                cactusWorkflowExperiment = CactusWorkflowExperiment(
                                                     sequences=self.sequences, 
                                                     newickTreeString=self.newickTree, 
                                                     #requiredSpecies=self.requiredSpecies,
                                                     #singleCopySpecies=self.singleCopySpecies,
                                                     databaseName="cactusAlignmentVanilla",
                                                     outputDir=tempLocalDir,
                                                     configFile=tempConfigFile)
                tempExperimentDir2 = os.path.join(tempLocalDir, "cactusAlignmentVanilla")
                cactusWorkflowExperiment.writeExperimentFile(tempExperimentFile2)
                
                
                #We're done with the progressive, now run the vanilla cactus for comparison
                tempJobTreeDir2 = os.path.join(tempLocalDir, "jobTreeVanilla")
                runCactusWorkflow(tempExperimentFile2, tempJobTreeDir2,
                                  jobTreeStats=True,
                                  setupAndBuildAlignments=True,
                                  buildReference=True,
                                  maxThreads=4)
                
                runJobTreeStatusAndFailIfNotComplete(tempJobTreeDir2)
                logger.info("Checked the job tree dir for the vanilla run")
                
                runCactusMAFGenerator(os.path.join(self.outputDir, "cactusVanilla.maf"), getCactusDiskString(tempExperimentDir2))
                
                system("jobTreeStats --jobTree %s --outputFile %s/jobTreeStatsVanilla.xml" % (tempJobTreeDir2, self.outputDir))
                system("mv %s %s" % (tempExperimentDir2, self.outputDir))
                system("mv %s %s/experimentVanilla.xml" % (tempExperimentFile2, self.outputDir))
            
            system("mv %s %s/config.xml" % (tempConfigFile, self.outputDir))
        #self.addChildTarget(MakeStats(outputFile, self.outputDir, self.options)) 

class MakeBlanchetteAlignments(Target):
    def __init__(self, options, params):
        Target.__init__(self)
        self.options = options
        self.params = params
        
    def run(self):
        outputDir = os.path.join(self.options.outputDir, "blanchette%s" % self.params)
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)
        repeats = 1
        for i in xrange(repeats):
            sequences, newickTreeString = getCactusInputs_blanchette(i)
            self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, os.path.join(outputDir, str(i)), 
                                        self.params))
        self.setFollowOnTarget(MakeBlanchetteStats(self.options, outputDir, repeats))
        
class MakeBlanchetteStats(Target):
    """Builds basic stats and the maf alignment.
    """
    def __init__(self, options, outputDir, repeats):
        Target.__init__(self)
        self.options = options
        self.outputDir = outputDir
        self.repeats = repeats
    
    def run(self):
        previousOutputFile = None
        previousOutputFile2 = None
        blanchettePath = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation")
        for i in xrange(self.repeats):
            trueAlignmentMFA = os.path.join(os.path.join(blanchettePath, "%.2i.job" % i), "true.mfa")
            trueAlignmentMAF = os.path.join(self.getLocalTempDir(), "temp.maf")
            treeFile = os.path.join(blanchettePath, "tree.newick")
            system("mfaToMaf --mfaFile %s --outputFile %s --treeFile %s" % (trueAlignmentMFA, trueAlignmentMAF, treeFile))
            trueRenamedMAF = os.path.join(trueAlignmentMAF, ".renamed")
            expPath = os.path.join(blanchettePath, "experiment.xml")
            applyNamingToMaf(expPath, trueAlignmentMAF, trueRenamedMAF)
            
            predictedAlignmentMaf = os.path.join(self.outputDir, str(i), "progressiveCactusAlignment", "Anc0", "Anc0.maf")
            
            outputFile = os.path.join(self.getLocalTempDir(), "temp%i" % i)
            system("mafComparator --mafFile1 %s --mafFile2 %s --outputFile %s" % (trueRenamedMAF, filteredPredictedAlignmentMAF, outputFile))
            system("cp %s %s" % (outputFile, os.path.join(self.outputDir, str(i), "mafComparison.xml")))
            if previousOutputFile != None:
                system("mergeMafComparatorResults.py --results1 %s --results2 %s --outputFile %s" % (outputFile, previousOutputFile, outputFile))
            previousOutputFile = outputFile
            
        system("mv %s %s" % (previousOutputFile, os.path.join(self.outputDir, "mafComparison.xml")))   
        
class MakeEvolverPrimatesLoci1(MakeBlanchetteAlignments):
    def setupStats(self, outputDir, simDir):
        #Setup the stat computation
        trueMaf = os.path.join(simDir, "all.burnin.maf")
        predictedMaf = os.path.join(outputDir, "progressiveCactusAlignment", "Anc0", "Anc0.maf")
        outputFile = os.path.join(outputDir, "mafComparison.burnin.xml")
        self.setFollowOnTarget(MakeStats(self.options, trueMaf, predictedMaf, outputFile))
    
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "loci1")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.chr6", "simChimp.chr6", "simGorilla.chr6", "simOrang.chr6"))
        outputDir = os.path.join(self.options.outputDir, "evolverPrimatesLoci1%s"  % self.params)
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, simDir)
        
        
class MakeEvolverMammalsLoci1(MakeEvolverPrimatesLoci1):
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "loci1")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.chr6", "simMouse.chr6", "simRat.chr6", "simCow.chr6", "simDog.chr6"))
        outputDir = os.path.join(self.options.outputDir, "evolverMammalsLoci1%s"  % self.params)
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, simDir)
        
class MakeStats(Target):
    def __init__(self, options, trueMaf, predictedMaf, outputFile):
        Target.__init__(self)
        self.options = options
        self.trueMaf = trueMaf
        self.predictedMaf = predictedMaf
        self.outputFile = outputFile
    
    def run(self):
        if not os.path.exists(self.outputFile):
            outputFile = os.path.join(self.getLocalTempDir(), "temp.xml")
            system("mafComparator --mafFile1 %s --mafFile2 %s --outputFile %s" % (self.trueMaf, self.predictedMaf, outputFile))
            system("mv %s %s" % (outputFile, self.outputFile))
      
class MakeAllAlignments(Target):
    """Makes alignments using pipeline.
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options
    
    def run(self):
        for params in SmallProgressive().generate():
            self.addChildTarget(MakeBlanchetteAlignments(self.options, params))
            #self.addChildTarget(MakeEvolverPrimatesLoci1(self.options, params))
            #self.addChildTarget(MakeEvolverMammalsLoci1(self.options, params))
            
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
    from progressiveBenchmarks.src.pipeline import *
    _test()
    main()

