######################################
#
#SequenceFinder finds all possible sequences of a given length (default 300) that contain SNPs,
#calculates the average mutual information value for all the SNPs in each sequence,
#sorts the sequences from lowest to highest avg mutual information value, and then 
#culls the sequences based on that value so that only the %50th percentile is kept.
#
#The sequences assume there's a primer at each end (default 15), so sequences only 
#look for SNPs within an area equal to the sequence length - 2 * primer size
#
#Data stored within each sequence is:
#   startPosition - the location of the first nucleotide in the sequence
#   sequenceLength - the total length of the sequence
#   avgMutualInformation - the mean mutual information, when comparing each SNP in
#           the sequence to every other SNP (sequences with 1 SNP are automatically 
#           given a value of 1)
#   SNPList - a list of SNP objects for the SNPs that fall within the sequence
#
######################################

from ISGDataPuller import SNP
import EntropyCalculator

# main function finds all sequences in ISGData of length sequenceLength and then culls those sequences based on their mutual information
def main(TreeTable, ISGData, sequenceLength = 300, primerSize = 15):

    sequences = []
    for i, SNP in enumerate(ISGData):
        sequence = Sequence(ISGData, i, sequenceLength, primerSize, TreeTable)
        sequenceIndex = 0
        if len(sequences):
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
