outputDir = output
maxCpus=3 
maxCpusPerAlignment=3
jobTreeParameters = --logDebug --maxThreads=${maxCpus} --maxJobs=${maxCpus}
batchSystemForAlignments = singleMachine singleMachine 1000000
databaseHost = localhost
databaseDir = dummy
parasolCommandForAlignment = parasol

all :
	#Nothing currently

run : 
	rm -rf jobTree
	python src/pipeline.py --parasolCommandForAlignment ${parasolCommandForAlignment} --databaseDir ${databaseDir} --batchSystemForAlignments "${batchSystemForAlignments}" --databaseHost ${databaseHost} --outputDir ${outputDir} --jobTree jobTree ${jobTreeParameters} --cpus=${maxCpusPerAlignment}
	rm -rf jobTree

clean :
	rm -rf ${outputDir}/*
	
test :
	#Nothing currently
