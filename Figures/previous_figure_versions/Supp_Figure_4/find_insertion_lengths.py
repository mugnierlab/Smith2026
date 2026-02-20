import statistics

reads = {}

with open("/Users/jaclyn.smith/labNotebook/erin_data/" + "VSG8_cleaned_recombinations.tsv", "r") as input_file:
    for line in input_file:
        if line.startswith("original"):
            pass
        else:
            try:
                reads[line.split("\t")[0]].append(line.strip())
            except KeyError:
                reads[line.split("\t")[0]] = [line.strip()]

insertion_list = {}
for entry in reads.keys():
    if len(reads[entry])%2 == 0:
        # even number of cross overs
        if reads[entry][0].split("\t")[6] == "target_to_donor":
            start = "donor"
            start_pos = int(reads[entry][0].split("\t")[5])
        else:
            # this doesn't occur
            pass
        i = 1
        for insertion in range(0, len(reads[entry]), 2):
            insert_start = reads[entry][insertion].split("\t")[5]
            insert_end = reads[entry][insertion + 1].split("\t")[4]
            insert_len = int(insert_end) - int(insert_start)
            insertion_list[entry + "_" + str(i)] = "\t".join(reads[entry][0].split("\t")[0:3]) + "\t" + str(insert_len)
            i += 1
        pass
    else:
        # odd number of cross overs
        if reads[entry][0].split("\t")[6] == "target_to_donor":
            # this doesn't occur
            pass
        else:
            insertion_list[entry + "_" + str(1)] = "\t".join(reads[entry][0].split("\t")[0:3]) + "\t" + reads[entry][0].split("\t")[4]
            i = 2
            for insertion in range(1, len(reads[entry]), 2):
                insert_start = reads[entry][insertion].split("\t")[5]
                insert_end = reads[entry][insertion + 1].split("\t")[4]
                insert_len = int(insert_end) - int(insert_start)
                insertion_list[entry + "_" + str(i)] = "\t".join(reads[entry][0].split("\t")[0:3]) + "\t" + str(insert_len)
                i += 1

with open("insert_len_VSG8_nanopore.tsv", "w") as output_file:
    insert_len_list = []
    for entry in insertion_list:
        output_file.write(entry + "\t" + insertion_list[entry] + "\n")
        insert_len_list.append(int(insertion_list[entry].split("\t")[-1]))

print(statistics.mean(insert_len_list))
print(statistics.median(insert_len_list))
