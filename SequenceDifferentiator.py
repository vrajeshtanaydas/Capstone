def main():
    print("tryopen")
    input_file = "./TestISG.txt"
    fopen = open(input_file, "r")
    print("opned")
    sequence_groups = []
    sequence_length = 300
    
    # skip genome number
    fopen.seek(12)
    num_genomes = int(fopen.readline().strip())
    num_genomes += 1

    # read header line with species names
    species_names = fopen.readline().strip().split('\t')
    name_counter = 0
    for name in species_names:
        if name_counter > 1 and name_counter < (num_genomes + 2):
            sequence_groups.append(name)
        name_counter += 1

    # get the sequence window data
    pos_counter = 0
    pos_start = None

    # while restriction
    while True:
        group_data = []
        break_while = False
        print("BLSJDLBJ")
        pos = fopen.readline().strip().split('\t')
        pos_start = int(pos[1])
        pos_end = pos_start + sequence_length
    
        # get SNP data (each line)
        for line in fopen.readlines():

            species_data = line.strip().split('\t')
        
            if int(species_data[1]) > pos_end or int(species_data[1]) < pos_start:
                break_while = True
                break

            group_data.append(species_data[1])
            group_data.append("")
            for snp_counter in range(2, num_genomes+2):
                group_data.append(species_data[snp_counter])

            pos_counter += 1

        if break_while is False:
            print("Break")
            break

        # do calculations
        matched_genomes = []
        for i in range(0, len(group_data)-1):
            i_matches = []
            if group_data[i] in matched_genomes:
                continue
            else:
                matched_genomes.append(group_data[i])
                i_matches.append(group_data[i])

            for j in range(i+1, len(group_data)):
                if group_data[i] == group_data[j]:
                    matched_genomes.append(group_data[j])
                    i_matches.append(group_data[j])

        print("BLHSDFI") 
        out_file = open("./differentiator_output.txt", "w+")
        for match in matched_genomes:
            out_file.write(match + "\n")


main()
