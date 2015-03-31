from ISGDataPuller import SNP
import EntropyCalculator

# main function finds all sequences in ISGData of length sequenceLength and then culls those sequences based on their mutual information
def main(TreeTable, ISGData, sequenceLength = 300, primerSize = 15):

    sequences = []
    for i, SNP in enumerate(ISGData):
        sequence = Sequence(ISGData, i, sequenceLength, primerSize, TreeTable)
        sequenceIndex = 0
        if len(sequences):
            #print(len(sequences))
            while sequence.avgMutualInformation >= sequences[sequenceIndex].avgMutualInformation:
                sequenceIndex += 1
                if sequenceIndex == len(sequences):
                    break
        sequences.insert(sequenceIndex, sequence)
    
    #culling sequences, might need to rework 
    mutualInformationRange = sequences[-1].avgMutualInformation - sequences[0].avgMutualInformation
    mutualInformationThreshold = mutualInformationRange / 2 + sequences[0].avgMutualInformation
    thresholdIndex = 0
    while sequences[thresholdIndex].avgMutualInformation <= mutualInformationThreshold and thresholdIndex < len(sequences) - 1:
        thresholdIndex += 1
    sequences = sequences[:thresholdIndex]

    return sequences

# class containing information about a specific sequence
class Sequence(object):


    def __init__(self, ISGData, snpIndex, sequenceLength, primerSize, TreeTable):
        self.avgMutualInformation = 1
        self.sequenceLength = sequenceLength
        self.startPosition = ISGData[snpIndex].pos - primerSize
        if self.startPosition < 0:
            self.startPosition = 1
        self.SNPList = []
        position = self.startPosition
        # adds SNPs to sequence as long as their position is within the sequenceLength
        while position <= (self.startPosition + sequenceLength - primerSize):

            self.SNPList.append(ISGData[snpIndex])
            snpIndex += 1
            if snpIndex == len(ISGData):
                break 
            position = ISGData[snpIndex].pos
        if len(self.SNPList) > 1:
            self.calcMutualInformation(TreeTable)
        else:
            print(len(self.SNPList))

    
    # compares each pair of SNPs to determine mutual information
    def calcMutualInformation(self, TreeTable):
        sumMutualInformation = 0
        pairCount = 0

        for i, snp in enumerate(self.SNPList):
            if i < len(self.SNPList) - 1:
                for j in range(i + 1, len(self.SNPList)):
                    hx = snp.entropy
                    hy = self.SNPList[j].entropy
                    hxy = EntropyCalculator.main(TreeTable, snp.aGenomes,snp.tGenomes,snp.cGenomes,snp.gGenomes, self.SNPList[j].aGenomes, self.SNPList[j].tGenomes,  self.SNPList[j].cGenomes,  self.SNPList[j].gGenomes)
                    sumMutualInformation += hx + hy - hxy
                    pairCount += 1

        self.avgMutualInformation = sumMutualInformation / pairCount
        print(self.avgMutualInformation)
