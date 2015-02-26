 
# since this is part of the pipeline, we need sys module
import sys
import EntropyCalculator

def main(inputFile, treeTable):

    # TODO -- need to handle exceptions better:
    #   file I/O?
    #   wrong file type input?

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

    def getChrom(self):
        return self.chrom

    def getPos(self):
        return self.pos

    def getEntropy(self):
        return self.entropy

    def getAGenomes(self):
        return self.aGenomes

    def getTGenomes(self):
        return self.tGenomes

    def getCGenomes(self):
        return self.cGenomes

    def getGGenomes(self):
        return self.gGenomes

    def printSNP(self):
        print(self.pos, self.entropy, self.aGenomes, self.tGenomes, self.cGenomes, self.gGenomes)



