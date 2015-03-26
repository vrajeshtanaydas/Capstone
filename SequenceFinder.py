from ISGDataPuller import SNP

# main function finds all sequences in ISGData of length sequenceLength and then culls those sequences based on their mutual information
def main(ISGData, sequenceLength = 300, primerSize = 15):
    sequences = []
    for i, SNP in enumerate(ISGData):
        sequence = Sequence(ISGData, i, sequenceLength - 2*primerSize)
        sequenceIndex = 0
        if len(sequences):
            print(len(sequences))
            while sequence.getavgMutualInformation() >= sequences[sequenceIndex].getavgMutualInformation():
                sequenceIndex += 1
                if sequenceIndex == len(sequences):
                    break
        sequences.insert(sequenceIndex, sequence)
    
    #culling sequences, might need to rework 
    mutualInformationRange = sequences[-1].getavgMutualInformation() - sequences[0].getavgMutualInformation()
    mutualInformationThreshold = mutualInformationRange / 2 + sequences[0].getavgMutualInformation()
    thresholdIndex = 0
    while sequences[thresholdIndex].getavgMutualInformation() <= mutualInformationThreshold and thresholdIndex < len(sequences) - 1:
        thresholdIndex += 1
    sequences = sequences[:thresholdIndex]
    
    return sequences

# class containing information about a specific sequence
class Sequence(object):


    def __init__(self, ISGData, snpIndex, sequenceLength):
        self.avgMutualInformation = 1
        self.startPosition = ISGData[snpIndex].getPos() - primerSize
        if self.startPosition < 0:
            self.startPosition = 1
        self.SNPList = []
        position = self.startPosition
        # adds SNPs to sequence as long as their position is within the sequenceLength
        while position <= (self.startPosition + sequenceLength - primersize):
            self.SNPList.append(ISGData[snpIndex])
            snpIndex += 1
            if snpIndex == len(ISGData):
                break 
            position = ISGData[snpIndex].getPos()
        #self.calcMutualInformation()
    
    def getSNPList(self):
        return self.SNPList
    
    def getStartPosition(self):
        return self.startPosition
    
    def getavgMutualInformation(self):
        return self.avgMutualInformation
    
    # compares each pair of SNPs to determine mutual information
    def calcMutualInformation(self):
        sumMutualInformation = 0
        pairCount = 0
        for i, snp in enumerate(self.SNPList):
            if i < len(self.SNPList) - 1:
                for j in range(i + 1, len(self.SNPList)):
                    hx = snp.getEntropy()
                    hy = self.SNPList[j].getEntropy()
                    hxy = pairEntropyCalculator(snp, self.SNPList[j])
                    sumMutualInformation += hx + hy - hxy
                    pairCount += 1
        self.avgMutualInformation = sumMutualInformation / pairCount
    