import os

def main(sequences, output_file):
    outfile = open(output_file, "w")
    # get the genome names
    species = sequences[0].SNPList[0].aGenomes
    species += sequences[0].SNPList[0].tGenomes
    species += sequences[0].SNPList[0].cGenomes
    species += sequences[0].SNPList[0].gGenomes

    sequence_length = len(sequences)

    # create a dictionary with keys as species names
    genome_names = []
    for i in range(len(species)):
        genome_names.append(species[i])

    groups = dict.fromkeys(genome_names, "")

    # create vars for holding snp information
    a_groups = c_groups = t_groups = g_groups = []

    # iterate through each sequence region
    for sequence in sequences:    
        # get the snp information in the current region for each species
        for i in range(len(sequence.SNPList)):
            a_groups = sequence.SNPList[i].aGenomes
            c_groups = sequence.SNPList[i].cGenomes
            t_groups = sequence.SNPList[i].tGenomes
            g_groups = sequence.SNPList[i].gGenomes
            for a in a_groups:
                if i == 0:
                    groups[a] = 'A'
                else:
                    groups[a] += 'A'
            for c in c_groups:
                if i == 0:
                    groups[c] = 'C'
                else:
                    groups[c] += 'C'
            for g in g_groups:
                if i == 0:
                    groups[g] = 'G'
                else:
                    groups[g] += 'G'
            for t in t_groups:
                if i == 0:
                    groups[t] = 'T'
                else:
                    groups[t] += 'T'
        # calculate the groupings of the species according to snps
        diffs = [[""]]
        appends = False                   
        # iterate through each species name
        for name in genome_names:
            # get the snp information from the groups dictionary
            snps = groups[name]
            for i in range(len(diffs)):
                # add the snp and species name information to groupings
                if snps == diffs[i][0]:
                    diffs[i].append(name)
                    appends = True
            # get rid of empty element in groupings 
            if appends == False:    
                if diffs[0][0] == "":
                   diffs[0][0] = snps
                   diffs[0].append(name)
                else:
                    diffs.append([snps, name])

        # output information to file
        mutual_info = str(sequence.avgMutualInformation)
        start_pos = str(sequence.startPosition)
        end_pos = str(sequence.startPosition + sequence.sequenceLength)
        outfile.write("sequence region: " + start_pos + "-" + end_pos)
        outfile.write("; mutual information: " + mutual_info + '\n')
        for d in diffs:
            outfile.write("[ ")
            for i in d:
                outfile.write(i + " ") 
            outfile.write("]\n")
    outfile.close()
