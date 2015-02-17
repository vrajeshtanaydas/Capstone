# since this is part of the pipeline, we need sys module
import sys
import ISGDataPuller
import EntropyCalculator
from TreeTableMaker import TreeTable

# expect first argument for the input file
try:
    inputFile = sys.argv[1]
except Exception:
    raise Exception("No argument provided as input")

# expect second argument for the output file
try:
    outputFile = sys.argv[2]
except Exception:
    raise Exception("No argument provided as output")

try:
    treeFile = sys.argv[3]
except Exception:
    raise Exception("No argument provided for phylogenetic tree")

dataPullerOutput = "dataPullerOutput.txt"

treeTable = TreeTable(treeFile)

ISGDataPuller.main(inputFile, treeTable, dataPullerOutput)

#EntropyCalculator.main(treeTable, dataPullerOutput)
