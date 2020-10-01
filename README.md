IBR-Capstone
========

Increasing Biogeographical Resolution


INSTRUCTIONS

Run the following command: 

Python.exe Pipeline.py <inputFile>.txt <outputFile> <treeFile>.txt <referenceFile>.fasta

The pipeline requires four arguments:
	inputFile - the name/path of the ISG data file
	outputFile - the name/path used to create the output files
	treeFile - the name/path of the phylogenetic tree file
	referenceFile - the name/path of the FASTA file used as the 
		reference genome in the ISG data

The pipeline creates two files, <outputFile>.txt and <outputFile>.fasta, which will be saved into the currently open directory by default, or in the directory indicated by the outputFile path.
	<outputFile>.txt contains all the sequences and the genomes/groups 
		that they differentiate in a human-readable format.
	<outputFile>.fasta contains all the sequences in FASTA format.

