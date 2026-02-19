# env: plasmid_consol
import glob
import os

# convert from fastq to fasta
for file_name in glob.glob("*.fastq"):
    os.system("sed -n \'1~4s/^@/>/p;2~4p\' " + file_name + " > " +  file_name.split(".fastq")[0] + ".fa")

    # consolidate reads with cd-hit by 99%
    flags = "-c 0.99 -n 8 -M 16000 -T 16 -d 0"
    os.system("cd-hit-est -i " + file_name.split(".fastq")[0] + ".fa -o " +\
              file_name.split(".fastq")[0] + "_consol.fasta " + flags)

    # read in sequences
    read_dict = {}
    with open(file_name.split(".fastq")[0] + ".fa", "r") as input_fasta:
        for line in input_fasta:
            if line.startswith(">"):
                read_name = line.strip().split(">")[1]
            else:
                sequence = line.strip()
                read_dict[read_name] = sequence
    read_dict[read_name] = sequence

    clusters = {}
    first = True
    # get consol_groups
    with open(file_name.split(".fastq")[0] + "_consol.fasta.clstr", "r") as input_cluster_file:
        clusters = {}
        read_list = []
        for line in input_cluster_file:
            if "Cluster" in line:
                if not first:
                    clusters[leader_name] = read_list
                read_list = []
                first = False
            elif "at" not in line:
                leader_name = line.split(">")[1].split("...")[0]
                read_list.append(line.split(">")[1].split("...")[0])
            else:
                read_list.append(line.split(">")[1].split("...")[0])
        clusters[leader_name] = read_list

    # subset cluster files, only those with 3 or more sequences
    with open(file_name.split(".fastq")[0] + "_consensus.fa", "w") as output_con:
        for item in clusters.keys():
            if len(clusters[item]) <= 10:
                pass
            else:
                with open("tmp_leader.fasta", "w") as output_tmp_leader:
                    output_tmp_leader.write(">" + item + "\n")
                    output_tmp_leader.write(read_dict[item] + "\n")
                with open("tmp_reads.fasta", "w") as output_tmp_fasta:
                    for entry in clusters[item]:
                        output_tmp_fasta.write(">" + entry + "\n")
                        output_tmp_fasta.write(read_dict[entry] + "\n")

                os.system("minimap2 -a --sam-hit-only -t 16 tmp_leader.fasta tmp_reads.fasta > output.sam")
                os.system("samtools view -bS output.sam > output.bam")
                os.system("samtools sort output.bam -o output.sorted.bam")
                os.system("rm *.fai")
                os.system("samtools mpileup -f tmp_leader.fasta output.sorted.bam > output_pileup.txt")

                consensus_seq = ""
                num_skips = 0
                coverage_count = 0
                write = True
                with open("output_pileup.txt", "r") as mpileup:
                    for line in mpileup:
                        if num_skips > 0:
                            num_skips -= 1
                            continue

                        # count disagreements
                        if int(line.split("\t")[3]) == 1:
                            coverage_count += 1
                            if coverage_count >= 500:
                                # i believe this is too many disagreements among the reads
                                # move to the next sequence and do not make a consensus
                                write = False
                                break

                        # using error tolerance to determine base at each position.
                        if int(line.split("\t")[3]) * 0.85 <= line.split("\t")[4].count("*"):
                            continue

                        elif int(line.split("\t")[3]) * 0.85 >= line.split("\t")[4].count(".") + line.split("\t")[4].count(","):
                            new_char = max(set(line.split("\t")[4].upper()), key=line.split("\t")[4].upper().count)
                            if new_char.isalnum():
                                try:
                                    int(new_char)
                                    consensus_seq += line.split("\t")[2]
                                except ValueError:
                                    consensus_seq += new_char
                            else:
                                consensus_seq += line.split("\t")[2]

                        elif line.split("\t")[4].count("-") >= int(line.split("\t")[3]) * 0.55:
                            numeric_counts = [x for x in line.split("\t")[4] if x.isnumeric()]
                            final_numeric_counts = "".join(numeric_counts)
                            numeric_counts_max = max(final_numeric_counts, key=final_numeric_counts.count)
                            num_skips = int(numeric_counts_max)

                        elif line.split("\t")[4].count("+") >= int(line.split("\t")[3]) * 0.55:
                            numeric_counts = [x for x in line.split("\t")[4] if x.isnumeric()]
                            final_numeric_counts = "".join(numeric_counts)
                            numeric_counts_max = max(final_numeric_counts, key=final_numeric_counts.count)
                            extra_seq = line.split("\t")[4].split(numeric_counts_max)[1][0:int(numeric_counts_max)]
                            consensus_seq += extra_seq

                        else:
                            consensus_seq += line.split("\t")[2]

                if write:
                    output_con.write(">" + item + "\n")
                    output_con.write(consensus_seq + "\n")
                    print(consensus_seq)
