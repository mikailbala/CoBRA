# script to parse blast output files and filter based on percent identity, evalue, length, N50 contig along with coverage if necessary.

import statistics
from collections import Counter


def calc_n50(list_of_lengths):
    """Calculate N50 for a sequence of numbers.

    Args:
        list_of_lengths (list): List of numbers.

    Returns:
        float: N50 value.

    """
    tmp = []
    for tmp_number in set(list_of_lengths):
        tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()

    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]

    return median


#input file and lists for later use.
input_file = input("Please enter your blast out file: ")
# list for filtered node names, endo targeted.
nodes_names = []
# temporary list for token separation of string.
tokens = []
# list for filtered node headers themselves, endo targeted.
filtered_nodes = []
# list contig_length
contig_lens = []
# list for alpha filtered nodes.
alpha_nodes = []


# open the file and iterate through the lines to parse metadata.
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

            # append the lengths of each contig to a list before filter to calculate n50.
            contig_lens.append(node_len)

            # filter to create a list of bacterial nodes.
            if perc_idn > 97.5 and node_len > 1000.0 and e_value < 0.010:
                # print(node_name, node_len, e_value)
                # grab the lines with e-value < 10 and node length > 1000 bases.
                filtered_nodes.append(current_line)
                # grab those nodes in a list
                nodes_names.append(node_name)

        if not current_line:
            break

# sort the node names to make a unique list and print them.
unq_node_names = (list(set(nodes_names)))
print("len of unique node names list:")
print(len(unq_node_names))
print("Unique node Names:")
for i in unq_node_names:
    print(i)


# list of unique node lengths
unq_node_lens = []
# calculate sum of those unique contigs:
for node in unq_node_names:
    tmp_tokens = node.split("\t")
    tmp_tokens2 = tmp_tokens[0].split('_')
    node_len = int(tmp_tokens2[3])
    unq_node_lens.append(node_len)

unq_node_sum = sum(unq_node_lens)
print("\nsum of unique node lengths: ", unq_node_sum)


# calculate n50
contig_N50 = calc_n50(contig_lens)
print("\nN50:", contig_N50)


# make a list of coverages for coverage filter.
coverages = []
# loop for alpha filter.
for node in filtered_nodes:
    # separate the string by tab
    tmp_tokens = node.split("\t")
    # print(tmp_tokens)
    # break
    # separate the first token by underscore to extract metadata
    tmp_tokens2 = tmp_tokens[0].split('_')
    # print(tmp_tokens2)
    # break
    # name the tokens to extract metadata based on the location in header
    node_num = float(tmp_tokens2[1])
    depth_coverage = float(tmp_tokens2[5])
    node_len = int(tmp_tokens2[3])
    e_value = float(tmp_tokens[10])
    node_name = tmp_tokens[0]
    perc_idn = float(tmp_tokens[2])

    # append al the coverages for later use.
    coverages.append(depth_coverage)

    # this filter is specific to every library. please adjust according to the results.
    # if depth_coverage < 22.0: # manual filter - ccoL_pantoae_bc
    # if depth_coverage > 1.0: # manual filter - cbp
    if (node_len >= ((contig_N50 / 4) * 3)):  # and (depth_coverage < 4.0): # manual filter - ph_r_bc
        # if node_len >= (contig_N50 / 8) and depth_coverage < 3.0: # manual filter - ph_a_bc
        # if node_len >= contig_N50:
        alpha_nodes.append(node_name)

# make a unique list of double filtered nodes. and print them
unq_alpha = list(set(alpha_nodes))
print("\nlen of alpha node list:")
print(len(unq_alpha))
print("\nnodes:")
for i in unq_alpha:
    print(i)

# write to output file for circos.
out_node_file = input_file.replace("blastn.out", "alpha_nodes.txt")
with open(out_node_file, "w") as outfile:
    for i in unq_alpha:
        outfile.write(i)
        outfile.write("\n")


# list of alpha node lengths
unq_alpha_node_lens = []
# calculate sum of those alpha contigs:
for node in unq_alpha:
    tmp_tokens = node.split("\t")
    tmp_tokens2 = tmp_tokens[0].split('_')
    node_len = int(tmp_tokens2[3])
    unq_alpha_node_lens.append(node_len)

unq_alpha_sum = sum(unq_alpha_node_lens)
print("sum of unique alpha node lengths: ", unq_alpha_sum)

# # activate this if coverage filter is required.
# beta_nodes = []
#
#
# # loop through filtered_nodes to filter based on coverage. result: beta_nodes.
# for node in filtered_nodes:
#     # separate the string by tab
#     tmp_tokens = node.split("\t")
#     # separate the first token by underscore to extract metadata
#     tmp_tokens2 = tmp_tokens[0].split('_')
#
#     # name the tokens to extract metadata based on the location in header
#     node_num = float(tmp_tokens2[1])
#     depth_coverage = float(tmp_tokens2[5])
#     node_len = int(tmp_tokens2[3])
#     e_value = float(tmp_tokens[10])
#     node_name = tmp_tokens[0]
#     perc_idn = float(tmp_tokens[2])
#
#     # select a coverage of your choice.
#     if depth_coverage > 1.0:
#         beta_nodes.append(node_name)
#
# # make a unq list of beta nodes and print them.
# unq_beta = list(set(beta_nodes))
# print("\nlen of beta node list:")
# print(len(unq_beta))
# print("\nnodes:")
# for i in unq_beta:
#     print(i)
