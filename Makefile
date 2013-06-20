outputDir = output

maxCpusPerAlignment=15
batchSystemForAlignments = singleMachine
extraJobTreeArgs = 
databaseHost = localhost
databaseDir = ""
parasolCommandForAlignment = parasol

#maxCpusPerAlignment=25
#batchSystemForAlignments = parasol
#extraJobTreeArgs = --bigBatchSystem singleMachine --bigMemoryThreshold 2147483648 --bigMaxMemory 536870912000 --bigCpuThreshold 4 --bigMaxCpus 25
#databaseHost = kolossus-10
#databaseDir = /hive/users/benedict/progressiveBenchmarks/databases
#parasolCommandForAlignment = /hive/users/benedict/parasol -host=swarm-10

jobTreeParameters = --logDebug

all :
	#Nothing currently

run : 
	rm -rf jobTree
	python src/pipeline.py --parasolCommandForAlignment "${parasolCommandForAlignment}" ${extraJobTreeArgs} --defaultMemory 8589934593 --maxThreads 1 --databaseDir ${databaseDir} --batchSystemForAlignments "${batchSystemForAlignments}" --databaseHost ${databaseHost} --outputDir ${outputDir} --jobTree jobTree ${jobTreeParameters} --cpus=${maxCpusPerAlignment}
	rm -rf jobTree

clean :
	rm -rf ${outputDir}/*
	
test :
	#Nothing currently

