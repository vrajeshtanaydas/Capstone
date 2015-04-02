####################################################
#
#ISGDatapuller opens an ISG data file and converts the information
#into a list of SNPs. Each SNP contains:
#   chrom - the reference id for the reference genome
#   pos - the position within the genome that the SNP is at
#   entropy - the entropy value of the SNP, used to calculate mutual information
#   aGenomes - a list of all genomes that have the value A in the SNP
#   tGenomes - a list of all genomes that have the value Tin the SNP
#   cGenomes - a list of all genomes that have the value C in the SNP
#   gGenomes - a list of all genomes that have the value G in the SNP
#
#The entropy value for each SNP is calculated when its created. 
#
####################################################

import EntropyCalculator

def main(inputFile, treeTable):

    fo = open(inputFile, "r")

    ISGData = []

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

        # calculate entropy values
        entropy = EntropyCalculator.main(treeTable, frozenset(arrA), frozenset(arrT), frozenset(arrC), frozenset(arrG))

        ISGData.append(SNP(strChrom, int(strPos), entropy, arrA, arrT, arrC, arrG))

            
    fo.close()

    return ISGData

# class containing information about a specific SNP
class SNP(object):

    def __init__(self, chrom, pos, entropy, aGenomes, tGenomes, cGenomes, gGenomes):
        self.chrom = chrom
        self.pos = pos
        self.entropy = entropy
        self.aGenomes = aGenomes
        self.tGenomes = tGenomes
        self.cGenomes = cGenomes
        self.gGenomes = gGenomes

    def printSNP(self):
        print(self.pos, self.entropy, self.aGenomes, self.tGenomes, self.cGenomes, self.gGenomes)
