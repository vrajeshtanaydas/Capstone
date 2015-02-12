# since this is part of the pipeline, we need sys module
import sys

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

# TODO -- need to handle exceptions better:
#   file I/O?
#   wrong file type input?

fo = open(inputFile, "r")
output = open(outputFile, "w+")

#read in number of genomes from file
fo.seek(12)
numGenomes = int(fo.readline().strip())
numGenomes += 1            #add one to include reference genome 

#read in header row
arrLine = fo.readline().strip().split("\t")

#place genome names in array
arrGenomeName = []
intCounter = 0
for word in arrLine:
    if(intCounter > 1 and intCounter < (numGenomes + 2)):
        arrGenomeName.append(word)
    intCounter += 1
    
#Determine genome group SNP differentiates 
for line in fo.readlines():
    arrLine = line.strip().split("\t")
    strChrom = arrLine[0]
    strPos = arrLine[1]
    arrSNP = []
    #place all SNPs into array
    for i in range(2,numGenomes+2):
        arrSNP.append(arrLine[i])
    #sort Genomes into groups by SNP call
    arrA = []
    arrT = []
    arrC = []
    arrG = []
    intCounter = 0
    for j in arrSNP:
        if(arrSNP[intCounter] == 'A'):
            arrA.append(arrGenomeName[intCounter])
        if(arrSNP[intCounter] == 'T'):
            arrT.append(arrGenomeName[intCounter])
        if(arrSNP[intCounter] == 'C'):
            arrC.append(arrGenomeName[intCounter])
        if(arrSNP[intCounter] == 'G'):
            arrG.append(arrGenomeName[intCounter])
        intCounter += 1
    #Checking if there at least 2 groups of at least 2 genomes
    intGroups = 0
    if(len(arrA) >= 1):
        intGroups += 1
    if(len(arrT) >= 1):
        intGroups += 1
    if(len(arrC) >= 1):
        intGroups += 1
    if(len(arrG) >= 1):
        intGroups += 1
    #Write to file
    if(intGroups >= 1):
        output.write(strChrom + "\t" + strPos + "\n")
        #Only write groups with more than 1 element
        if(len(arrA) >= 1):
            output.write("[")
            for strGenome in arrA:
                output.write(strGenome)
                output.write("\t")
            output.write("] \n")
        if(len(arrT) >= 1):
            output.write("[")
            for strGenome in arrT:
                output.write(strGenome)
                output.write("\t")
            output.write("] \n")
        if(len(arrC) >= 1):
            output.write("[")
            for strGenome in arrC:
                output.write(strGenome)
                output.write("\t")
            output.write("] \n")
        if(len(arrG) >= 1):
            output.write("[")
            for strGenome in arrG:
                output.write(strGenome)
                output.write("\t")
            output.write("] \n")
            
        output.write("\n")
        
fo.close()
output.close()
