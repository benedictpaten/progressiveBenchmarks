outputDir = output
maxCpus=4
maxCpusPerAlignment=4
jobTreeParameters = --logDebug --maxThreads=${maxCpus} --maxJobs=${maxCpus}

all :
	#Nothing currently

run : 
	rm -rf jobTree
	python src/pipeline.py --outputDir ${outputDir} --jobTree jobTree ${jobTreeParameters} --cpus=${maxCpusPerAlignment}
	rm -rf jobTree

clean :
	rm -rf ${outputDir}/*
	
test :
	#Nothing currently
