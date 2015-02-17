# since this is part of the pipeline, we need sys module
import sys
import ISGDataPuller

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

ISGDataPuller.main(inputFile, outputFile)

