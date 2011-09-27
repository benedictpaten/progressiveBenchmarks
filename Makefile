outputDir = output
jobTreeParameters = --logDebug

all : 
	rm -rf jobTree
	python src/pipeline.py --outputDir ${outputDir} --jobTree jobTree ${jobTreeParameters}
	rm -rf jobTree
