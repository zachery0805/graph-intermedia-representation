import os
import time
import sys
from sets import Set

#files = ["100-10", "100-50", "1000-10", "1000-250", "10000-50", "10000-2500", "100000-200", "100000-25000", "1000000-250", "1000000-250000"]
files = ["100-10"]

os.system("export OMP_NUM_THREADS=" + sys.argv[1])

print "Running tests\n------------------\n"

for filename in files:
    inputFilename = "input/" + filename
    outputFilename = "output/" + filename + "-par.bic"
    print "running bic " + filename
    start = time.time()
    os.system("./bic " + inputFilename + "  >" + outputFilename)
    end = time.time()
    print "\t... " + str(round(end-start, 5)) + " seconds"

print "\nConverting Output\n------------------\n"
for filename in files:
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

print "\nRunning Validation\n------------------\n"
for filename in files:
	print "verifying " + filename
	correctFilename = "output/" + filename + ".bic"
	testFilename = "output/c_" +  filename + "-par.bic"

	correctComponents = Set()
	currSet = Set()
	with open(correctFilename, "r") as correctFile:
		for line in correctFile:
			if line == "\n":
				correctComponents.add(currSet)
				currSet = Set()
			else:
				currSet.add(line.rstrip())

	foundComponents = Set()
	currSet = Set()
	with open(testFilename, "r") as testFile:
		for line in testFile:
			if line == "\n":
				foundComponents.add(currSet)
				currSet = Set()
			else:
				currSet.add(line.rstrip())

	if correctComponents == foundComponents:
		print "\t... PASSED"
	else:
		print "\t... FAILED"