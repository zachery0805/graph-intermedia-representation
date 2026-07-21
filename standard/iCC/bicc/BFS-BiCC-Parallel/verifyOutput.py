import sys
from sets import Set

filename = sys.argv[1]

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