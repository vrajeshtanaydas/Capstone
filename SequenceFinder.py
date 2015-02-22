def main(dataPullerOutput, sequenceLength):
    sequences = []
    fo = open(dataPullerOutput, 'r')
    lines = fo.readlines()
    for i, line in enumerate(lines):
        lines[i] = line.split('\t')
    lineIndex = 0
    while lineIndex < len(lines):
        sequence = Sequence(lines, lineIndex, sequenceLength)
        sequenceIndex = 0
        while sequence.getavgMutualInformation() >= sequences[sequenceIndex].getavgMutualInformation():
            sequenceIndex += 1
        sequences.insert(sequenceIndex, sequence)
        lineIndex += 5
    mutualInformationRange = sequences[-1].getavgMutualInformation() - sequences[0].getavgMutualInformation()
    mutualInformationThreshold = mutualInformationRange / 2 + sequences[0].getavgMutualInformation()
    thresholdIndex = 0
    while sequences[thresholdIndex].getavgMutualInformation() <= mutualInformationThreshold:
        thresholdIndex += 1
    sequences = sequences[:thresholdIndex]
    return sequences

class Sequence(object):


    def __init__(self, lines, lineIndex, sequenceLength):
        self.avgMutualInformation = 1
        self.startPosition = int(lines[lineIndex][1])
        self.SNPList = []
        position = self.startPosition
        while position <= (self.startPosition + sequenceLength):
            entropy = int(lines[lineIndex][2])
            self.SNPList.append((position, entropy))
            lineIndex += 5
            if lineIndex >= len(lines):
                break
            position = int(lines[lineIndex][1])
        self.calcMutualInformation()
    
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
                    hx = snp[1]
                    hy = self.SNPList[j][1]
                    hxy = pairEntropyCalculator(snp[0], self.SNPList[j][0])
                    sumMutualInformation += hx + hy - hxy
                    pairCount += 1
        self.avgMutualInformation = sumMutualInformation / pairCount
    