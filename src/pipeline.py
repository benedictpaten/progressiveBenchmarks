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
from sonLib.nxnewick import NXNewick
from sonLib.nxtree import NXTree 

from cactus.shared.test import getCactusInputs_blanchette

from jobTree.test.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

from cactus.shared.common import runCactusProgressive
from cactus.shared.common import runCactusCreateMultiCactusProject
from cactus.shared.test import getInputs

from sonLib.bioio import TestStatus
from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.configWrapper import ConfigWrapper
from cactus.progressive.cactus_createMultiCactusProject import cleanEventTree
from cactus.progressive.ktserverLauncher import KtserverLauncher
from progressiveBenchmarks.src.params import Params
from progressiveBenchmarks.src.paramsGenerator import ParamsGenerator
from progressiveBenchmarks.src.paramsGenerator import EverythingButSelf
from progressiveBenchmarks.src.paramsGenerator import AllProgressive
from progressiveBenchmarks.src.paramsGenerator import BasicProgressive
from progressiveBenchmarks.src.paramsGenerator import SmallProgressive
from progressiveBenchmarks.src.paramsGenerator import SingleCase
from progressiveBenchmarks.src.paramsGenerator import KyotoTycoon
from progressiveBenchmarks.src.paramsGenerator import LastzTuning, RepeatMasking
from progressiveBenchmarks.src.applyNamingToMaf import applyNamingToMaf
from progressiveBenchmarks.src.summary import Summary

def getRootPathString():
    """
    function for finding external location
    """
    import progressiveBenchmarks.src.pipeline
    i = os.path.abspath(progressiveBenchmarks.src.pipeline.__file__)
    return os.path.split(os.path.split(i)[0])[0] #os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

def getCactusWorkflowPathString():
    """
    function for finding external location of cactus workflow
    """
    import cactus.pipeline.cactus_workflow
    return os.path.split(os.path.abspath(cactus.pipeline.cactus_workflow.__file__))[0]
    
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
        Target.__init__(self)
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
        
    def getName(self):
        tok1, tok2 = os.path.split(self.outputDir)
        while '_' not in tok2:
            tok1, tok2 = os.path.split(tok1)
        return tok2
    
    def runProgressive(self):
        logger.debug("Going to put the alignment in %s" % self.outputDir)
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)
        
        if not os.path.exists(os.path.join(self.outputDir, "progressiveCactusAlignment")):
            xmlTree = ET.parse(os.path.join(getRootPathString(), "lib", "cactus_workflow_config.xml"))
            
            #Set the parameters
            tempLocalDir = os.path.join(self.outputDir, "tempProgressiveCactusAlignment")
            system("rm -rf %s" % tempLocalDir)
            os.mkdir(tempLocalDir)
            
            #Set the config parameters
            self.params.applyToXml(xmlTree)
            config = xmlTree.getroot()
            assert config is not None
            
            #Write the config file
            tempConfigFile = os.path.join(tempLocalDir, "config.xml")
            fileHandle = open(tempConfigFile, 'w')
            assert fileHandle is not None
            tree = ET.ElementTree(config)
            tree.write(fileHandle)
            fileHandle.close()
         
            #Make the experiment file
            tempExperimentFile = os.path.join(tempLocalDir, "experiment.xml")
            
            if self.options.databaseDir is not None and self.options.databaseDir != "":
                dbDir = self.options.databaseDir
                dbName = "temp_%s" % self.getName()
                system("rm -rf %s" % os.path.join(dbDir, dbName))
                if not os.path.isdir(dbDir):
                    os.mkdir(dbDir)
            else:
                dbDir = tempLocalDir
                dbName = "progressiveCactusAlignment"
                            
            if self.params.kyotoTycoon != False and self.params.kyotoTycoon != None:
                dbConfElem = ET.Element("st_kv_database_conf", type="kyoto_tycoon")
                if self.params.kyotoTycoon == True or str(self.params.kyotoTycoon).lower() == "true":
                    ktElem = ET.SubElement(dbConfElem, "kyoto_tycoon", host=self.options.databaseHost, port="1978", database_dir=dbDir)
                elif self.params.kyotoTycoon == "inMemory":
                    ktElem = ET.SubElement(dbConfElem, "kyoto_tycoon", host=self.options.databaseHost, port="1978", database_dir=dbDir, in_memory="true")
                elif self.params.kyotoTycoon == "inMemoryNoSnapshot":
                    ktElem = ET.SubElement(dbConfElem, "kyoto_tycoon", host=self.options.databaseHost, port="1978", database_dir=dbDir, in_memory="true", snapshot="false")                
            else:
                dbConfElem = None
            
            cactusWorkflowExperiment = CactusWorkflowExperiment(
                                                 sequences=self.sequences, 
                                                 newickTreeString=self.newickTree, 
                                                 #requiredSpecies=self.requiredSpecies,
                                                 #singleCopySpecies=self.singleCopySpecies,
                                                 databaseName=dbName,
                                                 outputDir=tempLocalDir,
                                                 configFile=tempConfigFile,
                                                 databaseConf = dbConfElem)
            cactusWorkflowExperiment.writeExperimentFile(tempExperimentFile)
            
            #The jobtree
            tempJobTreeDir = os.path.join(tempLocalDir, "jobTree")
            
            #The place to put the temporary experiment dir
            tempExperimentDir = os.path.join(tempLocalDir, "progressiveCactusAlignment")
            
      
            #The temporary experiment 
            runCactusCreateMultiCactusProject(tempExperimentFile, 
                                              tempExperimentDir)
            logger.info("Setup the cactus progressive experiment")
            
            #Make the logfile for jobTree
            jobTreeLogFile = os.path.join(self.outputDir, "log.txt")

            #hack to get the 20flys going
            event = "Anc0"
            if len(self.sequences) > 10:
                event = "Anc00"
            
            configWrapper = ConfigWrapper(config)
            joinMaf = configWrapper.getJoinMaf()
            
            runCactusProgressive(os.path.join(tempExperimentDir, "progressiveCactusAlignment_project.xml"), 
                                 tempJobTreeDir, 
                                 #batchSystem=batchSystem, 
                                 buildHal=True,
                                 #buildTrees=buildTrees, buildReference=buildReference,
                                 jobTreeStats=True,
                                 maxThreads=int(self.options.cpus),
                                 maxJobs=int(self.options.cpus),
                                 defaultMemory=int(self.options.defaultMemory),
                                 logLevel="CRITICAL",
                                 logFile = jobTreeLogFile,
                                 event=event,
                                 retryCount=2,
                                 batchSystem=self.options.batchSystemForAlignments,
                                 extraJobTreeArgumentsString="--parasolCommand '%s'" % self.options.parasolCommandForAlignment,
                                 profileFile=os.path.join(self.outputDir, "profileFile"))
            logger.info("Ran the progressive workflow")
            
            #Check if the jobtree completed sucessively.
            runJobTreeStatusAndFailIfNotComplete(tempJobTreeDir)
            logger.info("Checked the job tree dir for the progressive run")
            
            #Run the cactus tree stats
            if self.params.kyotoTycoon != "inMemoryNoSnapshot":
                expPath = os.path.join(tempExperimentDir, "Anc0", "Anc0_experiment.xml")
                exp = ExperimentWrapper(ET.parse(expPath).getroot())
                if exp.getDbType() == "kyoto_tycoon":
                    ktserver = KtserverLauncher()
                    ktserver.spawnServer(exp) 
                treeStatsFile = os.path.join(self.outputDir, "treeStats.xml")
                system("cactus_treeStats --cactusDisk \'%s\' --flowerName 0 --outputFile %s" %(exp.getDiskDatabaseString(),
                                                                                               treeStatsFile))
                if exp.getDbType() == "kyoto_tycoon":
                    ktserver.killServer(exp)
                
            #Now copy the true assembly back to the output
            system("mv %s %s/experiment.xml" % (tempExperimentFile, self.outputDir))
            system("mv %s %s" % (tempExperimentDir, self.outputDir))
            system("jobTreeStats --jobTree %s --outputFile %s/jobTreeStats.xml" % (tempJobTreeDir, self.outputDir))
            system("mv %s %s/config.xml" % (tempConfigFile, self.outputDir))
                
            #But keep a link to the multicactus project in its original path so we can navigate
            # the paths in the xml...
            actualResultsDir = os.path.join(os.path.abspath(self.outputDir), "progressiveCactusAlignment")
            tempResultsDir = os.path.join(self.outputDir, "tempProgressiveCactusAlignment")
            system("ln -s %s %s" % (actualResultsDir, tempResultsDir))
            
            #database dir given, so we erase it and overwrite with temp
            if self.options.databaseDir is not None and self.options.databaseDir != "":
                system("rm -rf %s && mv %s %s" % (os.path.join(dbDir, self.getName()),
                                                  os.path.join(dbDir, dbName),
                                                  os.path.join(dbDir, self.getName())))
                
    def runVanilla(self):
        logger.debug("Going to put the alignment in %s" % self.outputDir)
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)

        if not os.path.exists(os.path.join(self.outputDir, "cactusAlignmentVanilla")):
            xmlTree = ET.parse(os.path.join(getCactusWorkflowPathString(), "cactus_workflow_config.xml"))
            
            #Set the parameters
            tempLocalDir = os.path.join(self.outputDir, "tempVanillaCactusAlignment")
            system("rm -rf %s" % tempLocalDir)
            os.mkdir(tempLocalDir)
            
            #Set the config parameters
            #self.params.applyToXml(xmlTree)
            config = xmlTree.getroot()
            assert config is not None
        
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

            #Make the logfile for jobTree
            jobTreeLogFile = os.path.join(self.outputDir, "log.txt")
            
            cactusWorkflowExperiment = CactusWorkflowExperiment(
                                                 sequences=self.sequences, 
                                                 newickTreeString=self.newickTree, 
                                                 #requiredSpecies=self.requiredSpecies,
                                                 #singleCopySpecies=self.singleCopySpecies,
                                                 databaseName="cactusAlignmentVanilla",
                                                 outputDir=tempLocalDir,
                                                 configFile=tempConfigFile,
                                                 mafFile=os.path.join(self.outputDir, "cactusVanilla.maf"),
                                                 halFile=os.path.join(self.outputDir, "cactusVanilla.hal"),
                                                 fastaFile=os.path.join(self.outputDir, "cactusVanilla.fa"))
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
                              buildReference=True,
                              buildHal=True,
                              buildMaf=True,
                              buildFasta=True,
                              #Passing options
                              maxThreads=int(self.options.cpus),
                              maxJobs=int(self.options.cpus),
                              logFile=jobTreeLogFile)
            
            runJobTreeStatusAndFailIfNotComplete(tempJobTreeDir2)
            logger.info("Checked the job tree dir for the vanilla run")
            
            #runCactusMAFGenerator(os.path.join(self.outputDir, "cactusVanilla.maf"), getCactusDiskString(tempExperimentDir2))
            
            #Run the cactus tree stats
            treeStatsFile = os.path.join(self.outputDir, "treeStats.xml")
            system("cactus_treeStats --cactusDisk \'%s\' --flowerName 0 --outputFile %s" %(exp.getDiskDatabaseString(),
                                                                                        treeStatsFile))
            
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
            
            if self.params.vanilla == False:
                trueRenamedMAF = trueAlignmentMAF + ".renamed"
                expPath = os.path.join(self.outputDir, str(i), "experiment.xml")
                applyNamingToMaf(expPath, trueAlignmentMAF, trueRenamedMAF)
                trueAlignmentMAF = trueRenamedMAF            
                predictedAlignmentMaf = os.path.join(self.outputDir, str(i), "progressiveCactusAlignment", "Anc0", "Anc0.maf")
            else:
                predictedAlignmentMaf = os.path.join(self.outputDir, str(i), "cactusVanilla.maf")
            
            if os.path.exists(predictedAlignmentMaf):
                outputFile = os.path.join(self.getLocalTempDir(), "temp%i" % i)
                system("mafComparator --mafFile1 %s --mafFile2 %s --outputFile %s --sampleNumber 100000000 " % (trueAlignmentMAF, predictedAlignmentMaf, outputFile))
                system("cp %s %s" % (outputFile, os.path.join(self.outputDir, str(i), "mafComparison.xml")))
                if previousOutputFile != None:
                    system("mergeMafComparatorResults.py --results1 %s --results2 %s --outputFile %s" % (outputFile, previousOutputFile, outputFile))
                previousOutputFile = outputFile
        
        if previousOutputFile is not None and os.path.exists(previousOutputFile):    
            system("mv %s %s" % (previousOutputFile, os.path.join(self.outputDir, "mafComparison.xml")))   
        
class MakeEvolverPrimatesLoci1(MakeBlanchetteAlignments):
    name = "evolverPrimatesLoci1"
    def setupStats(self, outputDir, trueMaf, params):
        #Setup the stat computation
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
        self.setupStats(outputDir, os.path.join(simDir, "all.burnin.maf"), self.params)
        
        
class MakeEvolverMammalsLoci1(MakeEvolverPrimatesLoci1):
    name = "evolverMammalsLoci1"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "loci1")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.chr6", "simMouse.chr6", "simRat.chr6", "simCow.chr6", "simDog.chr6"))
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "all.burnin.maf"), self.params)
        
class MakeEvolverMammalsLoci1HumanMouse(MakeEvolverPrimatesLoci1):
    name = "evolverMammalsHumanMouseLoci1"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "loci1")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.chr6", "simMouse.chr6"))
        newickTreeString = "(simHuman:0.144018,simMouse:0.356483);"
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "all.burnin.maf"), self.params)
        
class MakeEvolverHumanMouseLarge(MakeEvolverPrimatesLoci1):
    name = "evolverHumanMouseLarge"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "large")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.trf.repmask.fa", "simMouse.trf.repmask.fa"))
        newickTreeString = "(simHuman:0.144018,simMouse:0.356483);"

        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "burnin.maf.map"), self.params)

class MakeHumanMouseWholeGenomes(MakeEvolverPrimatesLoci1):
    name = "humanMouseWholeGenomes"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "realMammals")
        sequences, newickTreeString = getInputs(simDir, ("hg19.fa.filterNs", "mm10.fa.filterNs"))
        newickTreeString = "(simHuman:0.144018,simMouse:0.356483);"
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "burnin.maf.map"), self.params)
        
class MakeHumanMouseDogWholeGenomes(MakeEvolverPrimatesLoci1):
    name = "humanMouseDogWholeGenomes"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "realMammals")
        sequences, newickTreeString = getInputs(simDir, ("hg19.fa.filterNs", "mm10.fa.filterNs", "canFam3.fa.filterNs"))
        newickTreeString = "((HUMAN:0.144018,MOUSE:0.356483)Anc0:0.0238,DOG:0.197)MRCA;"
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "burnin.maf.map"), self.params)

class MakeEvolverMammalsLarge(MakeEvolverPrimatesLoci1):
    name = "evolverMammalsLarge"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "large")
        sequences, newickTreeString = getInputs(simDir, ("simCow.trf.repmask.fa", "simDog.trf.repmask.fa", "simHuman.trf.repmask.fa", "simMouse.trf.repmask.fa", "simRat.trf.repmask.fa"))
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "burnin.maf.map"), self.params)
        
class MakeBlanchetteHumanMouse(MakeEvolverPrimatesLoci1):
    name = "blanchetteHumanMouse"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation", "00.job")
        sequences = os.path.join(simDir, "HUMANDIR"), os.path.join(simDir, "MOUSEDIR")
        #, newickTreeString = getInputs(simDir, ("HUMAN", "MOUSE"))
        newickTreeString = "(HUMAN:0.144018,MOUSE:0.356483);"
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "true.maf"), self.params)
        
class MakeEvolverMammalsLociMedium(MakeEvolverPrimatesLoci1):
    name = "evolverMammalsLociMedium"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "medium")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.fa", "simMouse.fa", "simRat.fa", "simCow.fa", "simDog.fa"))
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "burnin.maf"), self.params)
        
class MakeEvolverPrimatesMedium(MakeEvolverPrimatesLoci1):
    name = "evolverPrimatesMedium"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "medium")
        sequences, newickTreeString = getInputs(simDir, ("simGorilla.fa", "simHuman.fa", "simChimp.fa", "simOrang.fa"))
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "all.maf"), self.params)
   
class MakeEvolverPrimatesLarge(MakeEvolverPrimatesLoci1):
    name = "evolverPrimatesLarge"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "large")
        sequences, newickTreeString = getInputs(simDir, ("simGorilla.fa", "simHuman.fa", "simChimp.fa", "simOrang.fa"))
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "simPrimates.ancestor.maf"), self.params)


class Make20Flys(MakeEvolverPrimatesLoci1):
    name = "20Flys"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "flys", "20")
        sequences, newickTreeString = getInputs(simDir, ("droGri2.fa", "droVir3.fa", "droMoj3.fa", "droBip.fa", "droAna3.fa", "droKik.fa", "droFic.fa", "dm3.fa", "droSim1.fa", "droSec1.fa", "droYak2.fa", "droEre2.fa", "droEug.fa", "droBia.fa", "droTak.fa", "droEle.fa", "droRho.fa", "droPer1.fa", "dp4.fa", "droWil1.fa"))
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        #self.setupStats(outputDir, os.path.join(simDir, "burnin.maf.map"), self.params)

#########################
#####Tests added to for establishing correct repeat masking
#########################

class MakeBlanchetteHumanMouseDog(MakeEvolverPrimatesLoci1):
    name = "blanchetteHumanMouseDog"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation", "00.job")
        sequences = os.path.join(simDir, "HUMAN"), os.path.join(simDir, "MOUSE"), os.path.join(simDir, "DOG")
        #, newickTreeString = getInputs(simDir, ("HUMAN", "MOUSE"))
        newickTreeString = "((HUMAN:0.144018,MOUSE:0.356483)Anc0:0.0238,DOG:0.197)MRCA;"
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "true.maf"), self.params)

class MakeEvolverMammalsLoci1HumanMouseDog(MakeEvolverPrimatesLoci1):
    name = "evolverMammalsLoci1HumanMouseDog"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "loci1")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.chr6", "simMouse.chr6", "simDog.chr6"))
        newickTreeString = "((HUMAN:0.144018,MOUSE:0.356483)Anc0:0.0238,DOG:0.197)MRCA;" #Over-ride the full phylogeny
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "all.burnin.maf"), self.params)

class MakeEvolverMammalsLociMediumHumanMouseDog(MakeEvolverPrimatesLoci1):
    name = "evolverMammalsLociMediumHumanMouseDog"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "medium")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.fa", "simMouse.fa", "simDog.fa"))
        newickTreeString = "((HUMAN:0.144018,MOUSE:0.356483)Anc0:0.0238,DOG:0.197)MRCA;" #Over-ride the full phylogeny
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "burnin.maf"), self.params)

class MakeEvolverMammalsLargeHumanMouseDog(MakeEvolverPrimatesLoci1):
    name = "evolverMammalsLargeHumanMouseDog"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "large")
        sequences, newickTreeString = getInputs(simDir, ("simHuman.trf.repmask.fa", "simMouse.trf.repmask.fa", "simDog.trf.repmask.fa"))
        newickTreeString = "((HUMAN:0.144018,MOUSE:0.356483)Anc0:0.0238,DOG:0.197)MRCA;" #Over-ride the full phylogeny
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        self.setupStats(outputDir, os.path.join(simDir, "burnin.maf.map"), self.params)    
        
class Make3Worms(MakeEvolverPrimatesLoci1):
    name = "threeWorms"
    def run(self):
        simDir = os.path.join(TestStatus.getPathToDataSets(), "worms")
        sequences, newickTreeString = getInputs(simDir, ("cb4.fa", "caeRem4.fa", "caePb3.fa"))
        newickTreeString = "((cb4:0.462053,caeRem4:0.40839)Anc0:0.001,caePb3:0.5266235)MRCA;" # "((HUMAN:0.144018,MOUSE:0.356483)Anc0:0.0238,DOG:0.197)MRCA;" #Over-ride the full phylogeny
        outputDir = os.path.join(self.options.outputDir, "%s%s"  % (self.name, self.params))
        self.addChildTarget(MakeAlignment(self.options, sequences, newickTreeString, outputDir,
                                          self.params))
        #self.setupStats(outputDir, os.path.join(simDir, "burnin.maf.map"), self.params)

############End repeat masking tests
        
        
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
            if os.path.exists(self.predictedMaf):
                system("mafComparator --mafFile1 %s --mafFile2 %s --outputFile %s --sampleNumber 100000000" % (self.trueMaf, self.predictedMaf, outputFile))
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
        for testCategory in [MakeBlanchetteHumanMouse, MakeBlanchetteHumanMouseDog, MakeBlanchetteAlignments, 
                             MakeEvolverPrimatesLoci1, MakeEvolverMammalsLoci1HumanMouse, MakeEvolverMammalsLoci1,
                             MakeEvolverMammalsLociMedium, MakeEvolverPrimatesMedium, MakeEvolverHumanMouseLarge,
                             MakeEvolverMammalsLoci1HumanMouseDog, MakeEvolverMammalsLociMediumHumanMouseDog,
                             MakeEvolverMammalsLargeHumanMouseDog]:
            for name, i in self.getBaseNames(testCategory):
                summary = Summary()
                for params in self.paramsGenerator.generate():
                    rowName = name + str(params)
                    basePath = self.getPath(testCategory, params, i)
                    jobTreeStatsPath = os.path.join(basePath, "jobTreeStats.xml")
                    mafCompPath = os.path.join(basePath, "mafComparison.xml")
                    treeStatsPath = os.path.join(basePath, "treeStats.xml")
                    if params.vanilla is False:
                        projPath = os.path.join(basePath, "progressiveCactusAlignment", 
                                                "progressiveCactusAlignment_project.xml")
                    else:
                        projPath = None                    
                    if os.path.exists(mafCompPath):
                        summary.addRow(rowName, params, jobTreeStatsPath, mafCompPath, 
                                       treeStatsPath, projPath)
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
        #pg = EverythingButSelf()
        #pg = SingleCase()
        pg = KyotoTycoon()
        #pg = RepeatMasking()
        #pg = LastzTuning()
        for params in pg.generate():
            self.addChildTarget(MakeBlanchetteHumanMouse(self.options, params))
            #self.addChildTarget(MakeBlanchetteAlignments(self.options, params))
            #self.addChildTarget(MakeEvolverPrimatesLoci1(self.options, params))
            #self.addChildTarget(MakeEvolverMammalsLoci1HumanMouse(self.options, params))
            #self.addChildTarget(MakeEvolverMammalsLoci1(self.options, params))
            #self.addChildTarget(MakeEvolverMammalsLociMedium(self.options, params))
            #self.addChildTarget(MakeEvolverPrimatesMedium(self.options, params))
            #self.addChildTarget(MakeEvolverHumanMouseLarge(self.options, params))
            #self.addChildTarget(MakeEvolverPrimatesLarge(self.options, params))
            #self.addChildTarget(MakeEvolverMammalsLarge(self.options, params))
            #self.addChildTarget(Make20Flys(self.options, params))
            #self.addChildTarget(MakeHumanMouseWholeGenomes(self.options, params))
            #self.addChildTarget(MakeHumanMouseDogWholeGenomes(self.options, params))
            
            ###Repeat masking problems
            #self.addChildTarget(MakeBlanchetteHumanMouseDog(self.options, params))
            #self.addChildTarget(MakeEvolverMammalsLoci1HumanMouseDog(self.options, params))
            #self.addChildTarget(MakeEvolverMammalsLociMediumHumanMouseDog(self.options, params))
            #self.addChildTarget(MakeEvolverMammalsLargeHumanMouseDog(self.options, params))
            #self.addChildTarget(Make3Worms(self.options, params))
        
        self.setFollowOnTarget(MakeSummary(self.options, pg))

def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
    parser.add_option("--outputDir", dest="outputDir")
    parser.add_option("--cpus", dest="cpus")
    parser.add_option("--batchSystemForAlignments", default="singleMachine")
    parser.add_option("--databaseHost", default="localhost")
    parser.add_option("--databaseDir", default="")
    parser.add_option("--parasolCommandForAlignment", default="parasol")
    
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

