outputDir = output

maxCpusPerAlignment=30
batchSystemForAlignments = singleMachine
databaseHost = localhost
databaseDir = ""
parasolCommandForAlignment = parasol

#maxCpusPerAlignment=25
#batchSystemForAlignments = parasol singleMachine 2147483648 10000
#databaseHost = kolossus-10
#databaseDir = /hive/users/benedict/progressiveBenchmarks/databases
#parasolCommandForAlignment = /hive/users/benedict/parasol -host=swarm-10

jobTreeParameters = --logDebug

all :
	#Nothing currently

run : 
	rm -rf jobTree
	python src/pipeline.py --parasolCommandForAlignment "${parasolCommandForAlignment}" --defaultMemory 8589934593 --noCheckPoints --maxThreads 1 --databaseDir ${databaseDir} --batchSystemForAlignments "${batchSystemForAlignments}" --databaseHost ${databaseHost} --outputDir ${outputDir} --jobTree jobTree ${jobTreeParameters} --cpus=${maxCpusPerAlignment}
	rm -rf jobTree

clean :
	rm -rf ${outputDir}/*
	
test :
	#Nothing currently

