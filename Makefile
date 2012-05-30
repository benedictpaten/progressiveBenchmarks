outputDir = output
maxCpus=3 
maxCpusPerAlignment=3
batchSystemForAlignments = singleMachine
databaseHost = localhost
databaseDir = ""
parasolCommandForAlignment = parasol

#outputDir = output
#maxCpus=25 
#maxCpusPerAlignment=25
#batchSystemForAlignments = parasol singleMachine 2147483648 10000
#databaseHost = kolossus-10
#databaseDir = /hive/users/benedict/progressiveBenchmarks/databases
#parasolCommandForAlignment = /cluster/home/markd/pub/parasol -host=swarm-10

jobTreeParameters = --logDebug --maxThreads=${maxCpus} --maxJobs=${maxCpus}

all :
	#Nothing currently

run : 
	rm -rf jobTree
	python src/pipeline.py --parasolCommandForAlignment "${parasolCommandForAlignment}" --defaultMemory 4294967296 --databaseDir ${databaseDir} --batchSystemForAlignments "${batchSystemForAlignments}" --databaseHost ${databaseHost} --outputDir ${outputDir} --jobTree jobTree ${jobTreeParameters} --cpus=${maxCpusPerAlignment}
	rm -rf jobTree

clean :
	rm -rf ${outputDir}/*
	
test :
	#Nothing currently

