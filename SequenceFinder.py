def main(dataPullerOutput, sequenceLength):
    sequences = []
    fo = open(dataPullerOutput, 'r')
    lines = fo.readlines()
    for i, line in enumerate(lines):
        lines[i] = line.split('\t')
    lineIndex = 0
    while lineIndex < len(lines):
        sequences.append(Sequence(lines, lineIndex, sequenceLength))
        lineIndex += 5
    return sequences

class Sequence(object):


    def __init__(self, lines, lineIndex, sequenceLength):
        self.mutualInformation = 0
        self.startPosition = int(lines[lineIndex][1])
        self.SNPList = []
        position = self.startPosition
        while position <= (self.startPosition + sequenceLength):
            #entropy = int(lines[lineIndex][2])
            self.SNPList.append(position) # append (position, entropy)
            lineIndex += 5
            if lineIndex >= len(lines):
                break
            position = int(lines[lineIndex][1])
    
    def getSNPList(self):
        return self.SNPList
    
    def getStartPosition(self):
        return self.startPosition
    
    def setMutualInformation(self, value):
        self.mutualInformation = value
    
    def getMutualInformation(self):
        return self.mutualInformation