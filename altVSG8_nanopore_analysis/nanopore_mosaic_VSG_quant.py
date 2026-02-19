import glob
import os


def blast_against_VSGs(data_path):
    for file_name in glob.glob(data_path + "M*_2kb.fa"):
        updated_output_file = file_name.split(".fa")[0] + "_blastResults.tsv"
        os.system("blastn -query " + file_name + " -task \"blastn\" -db all_GC_VSGs.fasta -out " + \
                  updated_output_file + " -outfmt \"6 qseqid sseqid pident nident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore\" -num_threads 16 -max_target_seqs 1 -max_hsps 1")


def sort_VSG_hits(data_path):
    for file_name in glob.glob(data_path + "M*_2kb.fa"):
        updated_blast_file = file_name.split(".fa")[0] + "_blastResults.tsv"
        updated_nonVSG_file = file_name.split(".fa")[0] + "_nonVSGs.fa"
        updated_VSG_hits = file_name.split(".fa")[0] + "_VSGhit.fa"
        read_dict = {}
        with open(file_name, "r") as input_fasta:
            for line in input_fasta:
                if line.startswith(">"):
                    read_name = line.strip().split(">")[1].split()[0]
                else:
                    sequence = line.strip()
                    read_dict[read_name] = sequence
        with open(updated_blast_file, "r") as input_blast, open(updated_nonVSG_file, "w") as output_fasta, \
                open(updated_VSG_hits, "w") as output_vsg_fasta:
            for line in input_blast:
                if int(line.split("\t")[4]) <= 200:
                    output_fasta.write(">" + line.split("\t")[0] + "\n")
                    output_fasta.write(read_dict[line.split("\t")[0]] + "\n")
                else:
                    output_vsg_fasta.write(">" + line.split("\t")[0] + "\n")
                    output_vsg_fasta.write(read_dict[line.split("\t")[0]] + "\n")


def VSG_hit_counts(data_path):
    read_num = {}
    for item in glob.glob(data_path + "/M*_VSGhit.fa"):
        with open(item, "r") as input_file:
            count = 0
            for line in input_file:
                if line.startswith(">"):
                    count += 1
            read_num[item] = count

    with open(data_path + "VSG_count.tsv", "w") as output_file:
        for entry in read_num.keys():
            output_file.write("\t".join([entry, str(read_num[entry])]) + "\n")


def main():
    work_path = ""
    blast_against_VSGs(work_path)
    sort_VSG_hits(work_path)
    VSG_hit_counts(work_path)


if __name__ == '__main__':
    main()