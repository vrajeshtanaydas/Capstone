##############################################
#
#The pipeline requires four arguments:
#	inputFile - the name/path of the ISG data file
#	outputFile - the name used to create the output files
#	treeFile - the name/path of the phylogenetic tree file
#	referenceFile - the name/path of the FASTA file used as the 
#		reference genome in the ISG data
#
#The pipeline creates two files, <outputFile>.txt and <outputFile>.fasta.
#	<outputFile>.txt contains all the sequences and the genomes/groups 
#		that they differentiate in a human-readable format.
#	<outputFile>.fasta contains all the sequences in FASTA format.
#
#
##############################################


# since this is part of the pipeline, we need sys module
import sys
import ISGDataPuller
import EntropyCalculator
import SequenceFinder
import SequenceDifferentiator
from TreeMaker import TreeTable
import SequenceToFasta

# expect first argument for the input file
try:
    inputFile = sys.argv[1]
except Exception:
    raise Exception("No argument provided as input")
try:
    f = open(inputFile)
except Exception:
    raise Exception("Input file does not exist")
f.close()

# expect second argument for the output file
try:
    outputFile = sys.argv[2]
except Exception:
    raise Exception("No argument provided as output")

try:
    treeFile = sys.argv[3]
except Exception:
    raise Exception("No argument provided for phylogenetic tree")
try:
    f = open(treeFile)
except Exception:
    raise Exception("Tree file does not exist")
f.close()

try:
    referenceFile = sys.argv[4]
except Exception:
    raise Exception("No argument provided for reference file")
try:
    f = open(referenceFile)
except Exception:
    raise Exception("Reference file does not exist")
f.close()

treeTable = TreeTable(treeFile)

ISGData = ISGDataPuller.main(inputFile, treeTable)

sequences = SequenceFinder.main(treeTable,ISGData)

SequenceDifferentiator.main(sequences, outputFile)

SequenceToFasta.main(sequences, referenceFile, outputFile)
