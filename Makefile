outputDir = output
maxCpus=3 
maxCpusPerAlignment=3
jobTreeParameters = --logDebug --maxThreads=${maxCpus} --maxJobs=${maxCpus}
batchSystemForAlignments = singleMachine singleMachine 2147483648 1
databaseHost = localhost
databaseDir = ""
parasolCommandForAlignment = parasol

all :
	#Nothing currently

run : 
	rm -rf jobTree
	python src/pipeline.py --parasolCommandForAlignment "${parasolCommandForAlignment}" --databaseDir ${databaseDir} --batchSystemForAlignments "${batchSystemForAlignments}" --databaseHost ${databaseHost} --outputDir ${outputDir} --jobTree jobTree ${jobTreeParameters} --cpus=${maxCpusPerAlignment}
	rm -rf jobTree

clean :
	rm -rf ${outputDir}/*
	
test :
	#Nothing currently
