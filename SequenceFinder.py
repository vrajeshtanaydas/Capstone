from ISGDataPuller import SNP

def main(ISGData, sequenceLength):
    sequences = []
    for i, SNP in enumerate(ISGData):
        sequence = Sequence(ISGData, i, sequenceLength)
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
    while sequences[thresholdIndex].getavgMutualInformation() <= mutualInformationThreshold:
        thresholdIndex += 1
    sequences = sequences[:thresholdIndex]
    
    return sequences

class Sequence(object):


    def __init__(self, ISGData, snpIndex, sequenceLength):
        self.avgMutualInformation = 1
        self.startPosition = ISGData[snpIndex].getPos()
        self.SNPList = []
        position = self.startPosition
        while position <= (self.startPosition + sequenceLength):
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
    