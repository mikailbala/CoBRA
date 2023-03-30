# script to parse blast output files and filter based on percent identity, evalue, length, N50 contig along with coverage if necessary.

# this script is a filter test for pantoea carbhekki.

#input file and lists for later use.
# input_file = input("Please enter your blast out file: ")
input_file = "cbp_blastn.out"

# create the result file.
output_file = "PI_filter_test_results.tsv"

# print status info
print("initiating the PI threshold test...\n")

# list for unique node sums after every loop.
unq_node_sums = []
# list for percent recoveries.
perc_recoveries = [] 
# list for percent identities
perc_idens = []

# set a for loop for changing percent identity values.
for i in range(50,101):
    current_perc_idn = i

    # reset everything for every new loop
    # list for filtered node names, endo targeted.
    nodes_names = []
    # temporary list for token separation of string.
    tokens = []
    # list for filtered node headers themselves, endo targeted.
    filtered_nodes = []
    # list contig_length
    contig_lens = []
    # list of unique node lengths
    unq_node_lens = []
    # sum of unq node lens
    unq_node_sum = 0
    # reset percent recovery
    perc_rec = 0

    # print updated status info
    # print("current percent identity value is: ", current_perc_idn)
    
        # open the file and loop through the lines to parse metadata.
    with open(input_file, "r") as infile:
        while True:
            current_line = infile.readline()

            if len(current_line) > 0:
                # separate the string by tab
                tokens = current_line.split("\t")
                # separate the first token by underscore to extract metadata
                tokens2 = tokens[0].split('_')
                # print(tokens)
                # print(tokens2)
                # name the tokens to extract metadata based on the location in header
                node_num = float(tokens2[1])
                depth_coverage = float(tokens2[5])
                node_len = int(tokens2[3])
                e_value = float(tokens[10])
                node_name = tokens[0]
                perc_idn = float(tokens[2])
                al_len = float(tokens[3])

                # append the lengths of each contig to a list before filter to calculate n50.
                contig_lens.append(node_len)
                
                if (perc_idn > current_perc_idn) and \
                        (node_len > 1000.0) and \
                        (e_value < 0.010) and \
                        (al_len >= ((node_len / 4) * 3)):
                    # grab the lines that pass the filter.
                    filtered_nodes.append(current_line)
                    # grab those node names in a list to make a unique node name list later.
                    nodes_names.append(node_name)

            if not current_line:
                break


    # sort the node names to make a unique list and print them.
    unq_node_names = (list(set(nodes_names)))

    # calculate sum of those unique contigs:
    for node in unq_node_names:
        tmp_tokens = node.split("\t")
        tmp_tokens2 = tmp_tokens[0].split('_')
        node_len = int(tmp_tokens2[3])
        unq_node_lens.append(node_len)

    # calculate sum of node lens.
    unq_node_sum = sum(unq_node_lens)

    #set a default bacterial genome size depending on the test organism.
    # this case pantoea is 1.2mb
    bact_gnm_size = 1200000

    # calculate percent recovery.
    perc_rec = (unq_node_sum * 100)/bact_gnm_size

    # append all that parsed metadata to write in the output file.
    unq_node_sums.append(unq_node_sum)
    perc_recoveries.append(perc_rec)
    perc_idens.append(current_perc_idn)


print((perc_idens))
print((perc_recoveries))



# write the output file.
with open(output_file,"w") as outFile:
    for i in range(len(perc_idens)):
        outFile.write(str(perc_idens[i]) + "\t" + str(perc_recoveries[i]) + "\n")
