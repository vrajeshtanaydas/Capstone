#########################################
#
#SequenceToFasta converts each sequence into FASTA format in a text file,
#which can then be read by other programs. The output file is organized 
#into most valuable sequences, based on the average mutual information 
#value of the sequence. The reference file is the same reference file used 
#to generate the ISG data.
#
#########################################

from SequenceFinder import Sequence
from ISGDataPuller import SNP
from math import floor, ceil

def main(sequences, referenceFile):
    fo = open("sequences.fasta", "w")
    for sequence in sequences:
        # writes the ID and position
        fo.write(">" + sequence.SNPList[0].chrom + " Pos: " + str(sequence.startPosition) + "-" + str(sequence.startPosition + sequence.sequenceLength) + "\n")
        reference = open(referenceFile, "r")
        # skips the first line and stores its length
        offset = len(reference.readline())
        # reads the second line and stores its length without the new line
        lineLength = len(reference.readline()) - 1
        # checks for blank lines
        while not lineLength:
            lineLength = len(reference.readline()) - 1
        # seeks to the start of the sequence
        reference.seek(sequence.startPosition + floor(sequence.sequenceLength / lineLength) + offset - 1, 0)
        # reads in sequenceLength characters plus max new lines, removes the new lines, then cuts off extra characters
        sequenceString = reference.read(sequence.sequenceLength + ceil(sequence.sequenceLength / lineLength)).replace("\n", "")[:sequence.sequenceLength]
        reference.close()
        # prints sequence with appropriate new lines
        for i, char in enumerate(sequenceString):
            if not i % lineLength and i != 0:
                fo.write("\n")
            fo.write(char)
        fo.write("\n")
    fo.close()