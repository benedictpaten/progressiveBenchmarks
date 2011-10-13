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
from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.cactus_createMultiCactusProject import cleanEventTree
from progressiveBenchmarks.src.params import Params
from progressiveBenchmarks.src.paramsGenerator import ParamsGenerator
from progressiveBenchmarks.src.paramsGenerator import EverythingButSelf
from progressiveBenchmarks.src.paramsGenerator import AllProgressive
from progressiveBenchmarks.src.paramsGenerator import BasicProgressive
from progressiveBenchmarks.src.paramsGenerator import SmallProgressive
from progressiveBenchmarks.src.paramsGenerator import SingleCase
from progressiveBenchmarks.src.applyNamingToMaf import applyNamingToMaf
from progressiveBenchmarks.src.summary import Summary

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
    
    def runProgressive(self):
        logger.debug("Going to put the alignment in %s" % self.outputDir)
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)

        if not os.path.exists(os.path.join(self.outputDir, "progressiveCactusAlignment")):
            config = ET.parse(os.path.join(getRootPathString(), "lib", "cactus_workflow_config.xml")).getroot()
            assert config is not None
            
            #Set the parameters
            tempLocalDir = os.path.join(self.outputDir, "tempProgressiveCactusAlignment")
            system("rm -rf %s" % tempLocalDir)
            os.mkdir(tempLocalDir)
            
            #Set the config parameters
            self.params.applyToXml(config)
        
            #Write the config file
            tempConfigFile = os.path.join(tempLocalDir, "config.xml")
            fileHandle = open(tempConfigFile, 'w')
            assert fileHandle is not None
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
                                 maxThreads=4,
                                 logLevel="DEBUG")
            logger.info("Ran the progressive workflow")
            
            #Check if the jobtree completed sucessively.
            runJobTreeStatusAndFailIfNotComplete(tempJobTreeDir)
            logger.info("Checked the job tree dir for the progressive run")
            
            #Now copy the true assembly back to the output
            system("mv %s %s/experiment.xml" % (tempExperimentFile, self.outputDir))
            system("mv %s %s" % (tempExperimentDir, self.outputDir))
            system("jobTreeStats --jobTree %s --outputFile %s/jobTreeStats.xml" % (tempJobTreeDir, self.outputDir))
            system("mv %s %s/config.xml" % (tempConfigFile, self.outputDir))
            
            #But keep a link to the multicactus project in its original path so we can navigate
            # the paths in the xml...
            system("ln -s %s %s" % (os.path.abspath(self.outputDir), tempExperimentDir))
                
    def runVanilla(self):
        logger.debug("Going to put the alignment in %s" % self.outputDir)
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)

        if not os.path.exists(os.path.join(self.outputDir, "cactusAlignmentVanilla")):
            config = ET.parse(os.path.join(getRootPathString(), "lib", "cactus_workflow_config.xml")).getroot()
            assert config is not None
            
            #Set the parameters
            tempLocalDir = os.path.join(self.outputDir, "tempVanillaCactusAlignment")
            system("rm -rf %s" % tempLocalDir)
            os.mkdir(tempLocalDir)
            
            #Set the config parameters
            self.params.applyToXml(config)
        
            #Write the config file
            tempConfigFile = os.path.join(tempLocalDir, "config.xml")
            fileHandle = open(tempConfigFile, 'w')
            assert fileHandle is not None
            tree = ET.ElementTree(config)
            tree.write(fileHandle)
            fileHandle.close()
         
            #Make the experiment file
            tempExperimentFile = os.path.join(tempLocalDir, "experiment.xml")
            #Now do standard cactus..
            #Make the experiment file
            tempExperimentFile2 = os.path.join(tempLocalDir, "experiment.xml")

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
            
            # apply naming to the event tree to be consistent with progressive
            exp = ExperimentWrapper(ET.parse(tempExperimentFile2).getroot())
            cleanEventTree(exp)
            exp.writeXML(tempExperimentFile2)
            
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
            
            system("jobTreeStats --jobTree %s --outputFile %s/jobTreeStats.xml" % (tempJobTreeDir2, self.outputDir))
            system("mv %s %s" % (tempExperimentDir2, self.outputDir))
            system("mv %s %s/experiment.xml" % (tempExperimentFile2, self.outputDir))
        
    
    def run(self):
        if self.params.vanilla == True:
           self.runVanilla()
        else:
            self.runProgressive()

class MakeBlanchetteAlignments(Target):
    name = "blanchette"
    def __init__(self, options, params):
        Target.__init__(self)
        self.options = options
        self.params = params
        
    def run(self):
        outputDir = os.path.join(self.options.outputDir, "%s%s" % (self.name, self.params))
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)
    
        for i in xrange(self.options.blanchetteRepeats):
            sequences, newickTreeString = getCactusInputs_blanchette(i)
            self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, os.path.join(outputDir, str(i)), 
                                        self.params))
        self.setFollowOnTarget(MakeBlanchetteStats(self.options, outputDir, self.params))
        
class MakeBlanchetteStats(Target):
    """Builds basic stats and the maf alignment.
    """
    def __init__(self, options, outputDir, params):
        Target.__init__(self)
        self.options = options
        self.outputDir = outputDir
        self.params = params
    
    def run(self):
        previousOutputFile = None
        previousOutputFile2 = None
        blanchettePath = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation")
        for i in xrange(self.options.blanchetteRepeats):
            trueAlignmentMFA = os.path.join(os.path.join(blanchettePath, "%.2i.job" % i), "true.mfa")
            trueAlignmentMAF = os.path.join(self.getLocalTempDir(), "temp.maf")
            treeFile = os.path.join(blanchettePath, "tree.newick")
            system("mfaToMaf --mfaFile %s --outputFile %s --treeFile %s" % (trueAlignmentMFA, trueAlignmentMAF, treeFile))
            
            
            trueRenamedMAF = trueAlignmentMAF + ".renamed"
            expPath = os.path.join(self.outputDir, str(i), "experiment.xml")
            applyNamingToMaf(expPath, trueAlignmentMAF, trueRenamedMAF)
            trueAlignmentMAF = trueRenamedMAF
            if self.params.vanilla == False:            
                predictedAlignmentMaf = os.path.join(self.outputDir, str(i), "progressiveCactusAlignment", "Anc0", "Anc0.maf")
            else:
                predictedAlignmentMaf = os.path.join(self.outputDir, str(i), "cactusVanilla.maf")
            
            outputFile = os.path.join(self.getLocalTempDir(), "temp%i" % i)
            system("mafComparator --mafFile1 %s --mafFile2 %s --outputFile %s" % (trueAlignmentMAF, predictedAlignmentMaf, outputFile))
            system("cp %s %s" % (outputFile, os.path.join(self.outputDir, str(i), "mafComparison.xml")))
            if previousOutputFile != None:
                system("mergeMafComparatorResults.py --results1 %s --results2 %s --outputFile %s" % (outputFile, previousOutputFile, outputFile))
            previousOutputFile = outputFile
            
        system("mv %s %s" % (previousOutputFile, os.path.join(self.outputDir, "mafComparison.xml")))   
        
class MakeEvolverPrimatesLoci1(MakeBlanchetteAlignments):
    name = "evolverPrimatesLoci1"
    def setupStats(self, outputDir, simDir, params):
        #Setup the stat computation
        trueMaf = os.path.join(simDir, "all.burnin.maf")
        if self.params.vanilla == False:
            predictedMaf = os.path.join(outputDir, "progressiveCactusAlignment", "Anc0", "Anc0.maf")
        else:
            predictedMaf = os.path.join(outputDir, "cactusVanilla.maf")
        outputFile = os.path.join(outputDir, "mafComparison.xml")
        self.setFollowOnTarget(MakeStats(self.options, trueMaf, predictedMaf, outputFile, params))
    
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "loci1")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.chr6", "simChimp.chr6", "simGorilla.chr6", "simOrang.chr6"))
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, simDir, self.params)
        
        
class MakeEvolverMammalsLoci1(MakeEvolverPrimatesLoci1):
    name = "evolverMammalsLoci1"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "loci1")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.chr6", "simMouse.chr6", "simRat.chr6", "simCow.chr6", "simDog.chr6"))
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, simDir, self.params)
        
class MakeStats(Target):
    def __init__(self, options, trueMaf, predictedMaf, outputFile, params):
        Target.__init__(self)
        self.options = options
        self.trueMaf = trueMaf
        self.predictedMaf = predictedMaf
        self.outputFile = outputFile
        self.params = params
    
    def run(self):
        if not os.path.exists(self.outputFile):
            outputFile = os.path.join(self.getLocalTempDir(), "temp.xml")
           
            trueRenamedMAF = os.path.join(self.getLocalTempDir(), "true_renamed.maf") 
            outputDir = os.path.split(self.outputFile)[0]
            expPath = os.path.join(outputDir, "experiment.xml")
            applyNamingToMaf(expPath, self.trueMaf, trueRenamedMAF)
            self.trueMaf = trueRenamedMAF
            system("mafComparator --mafFile1 %s --mafFile2 %s --outputFile %s" % (self.trueMaf, self.predictedMaf, outputFile))
            system("mv %s %s" % (outputFile, self.outputFile))

class MakeSummary(Target):
    def __init__(self, options, paramsGenerator):
        Target.__init__(self)
        self.options = options
        self.paramsGenerator = paramsGenerator
    
    def getBaseNames(self, testCategory):
        if testCategory == MakeBlanchetteAlignments:
            for i in xrange(self.options.blanchetteRepeats):
                yield (testCategory.name + str(i), i)
        else:
            yield (testCategory.name, None)
    
    def getPath(self, testCategory, params, i):
        if i is not None:
            return os.path.join(self.options.outputDir, testCategory.name + str(params), str(i))
        else:
            return os.path.join(self.options.outputDir, testCategory.name + str(params))
     
    def run(self):
        for testCategory in [MakeBlanchetteAlignments, MakeEvolverPrimatesLoci1, MakeEvolverMammalsLoci1]:
            for name, i in self.getBaseNames(testCategory):
                summary = Summary()
                for params in self.paramsGenerator.generate():
                    rowName = name + str(params)
                    basePath = self.getPath(testCategory, params, i)
                    jobTreeStatsPath = os.path.join(basePath, "jobTreeStats.xml")
                    mafCompPath = os.path.join(basePath, "mafComparison.xml")
                    if params.vanilla is False:
                        projPath = os.path.join(basePath, "progressiveCactusAlignment", 
                                                "progressiveCactusAlignment_project.xml")
                    else:
                        projPath = None                    
                    summary.addRow(rowName, params, jobTreeStatsPath, mafCompPath, projPath)
                summary.write(os.path.join(self.options.outputDir, "%s_summary.csv" % name))
                   
class MakeAllAlignments(Target):
    """Makes alignments using pipeline.
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options
    
    def run(self):
        #pg = ParamsGenerator()
        #pg = BasicProgressive()
        #pg = AllProgressive()
        pg = EverythingButSelf()
        #pg = SingleCase()
        for params in pg.generate():
            self.addChildTarget(MakeBlanchetteAlignments(self.options, params))
            self.addChildTarget(MakeEvolverPrimatesLoci1(self.options, params))
            self.addChildTarget(MakeEvolverMammalsLoci1(self.options, params))
        
        self.setFollowOnTarget(MakeSummary(self.options, pg))

        
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
    parser.add_option("--outputDir", dest="outputDir")
    
    Stack.addJobTreeOptions(parser)
    
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    options.blanchetteRepeats = 1
    
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

