outputDir = output
jobTreeParameters = --logDebug --maxThreads=1

all :
	#Nothing currently

run : 
	rm -rf jobTree
	python src/pipeline.py --outputDir ${outputDir} --jobTree jobTree ${jobTreeParameters}
	rm -rf jobTree

clean :
	rm -rf ${outputDir}/*
	
test :
	#Nothing currently
