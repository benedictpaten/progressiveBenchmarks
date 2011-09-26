outputDir = output

all : 
	python src/pipeline.py --outputDir ${outputDir}
