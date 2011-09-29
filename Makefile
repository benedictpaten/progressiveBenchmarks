outputDir = output
jobTreeParameters = --logDebug --maxThreads=1

all : 
	rm -rf jobTree
	python src/pipeline.py --outputDir ${outputDir} --jobTree jobTree ${jobTreeParameters}
	rm -rf jobTree

clean :
	rm -rf ${outputDir}/*
