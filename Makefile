outputDir = output

maxCpusPerAlignment=25
batchSystemForAlignments = singleMachine
extraJobTreeArgs = --maxMemory 746870912000 --maxCpus 25 --maxThreads 25
databaseHost = localhost
databaseDir = ""
parasolCommandForAlignment = parasol

#maxCpusPerAlignment=25
#batchSystemForAlignments = parasol
#extraJobTreeArgs = --bigBatchSystem singleMachine --bigMemoryThreshold 8589934592 --bigMaxMemory 746870912000 --bigCpuThreshold 4 --bigMaxCpus 25 --maxThreads 25
#databaseHost = juggernaut-101
#databaseDir = ""
#parasolCommandForAlignment = /hive/users/benedict/parasol -host=ku.local

jobTreeParameters = --logDebug

all :
	#Nothing currently

run : 
	rm -rf jobTree
	python src/pipeline.py --parasolCommandForAlignment "${parasolCommandForAlignment}" --extraJobTreeArgs '${extraJobTreeArgs}' --defaultMemory 8589934593 --maxThreads 1 --databaseDir ${databaseDir} --batchSystemForAlignments "${batchSystemForAlignments}" --databaseHost ${databaseHost} --outputDir ${outputDir} --jobTree jobTree ${jobTreeParameters} --cpus=${maxCpusPerAlignment}
	rm -rf jobTree

clean :
	rm -rf ${outputDir}/*
	
test :
	#Nothing currently

