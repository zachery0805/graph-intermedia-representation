import sys

filename = sys.argv[1]

print "converting " + filename
inputFilename = "output/" + filename + "-par.bic"
outputFilename = "output/c_" + filename + "-par.bic"

components = {}
with open(inputFilename, "r") as inputFile:
	for line in inputFile:
		line = line.rstrip()
		line = line.split(" = ")
		label = line[1]
		if not label in components:
			components[label] = []
		
		components[label].append(line[0])

with open(outputFilename, "w") as outputFile:
	for value in components.itervalues():
		for edge in value:
			outputFile.write(edge + "\n")
		outputFile.write("\n")

print "\t... done"