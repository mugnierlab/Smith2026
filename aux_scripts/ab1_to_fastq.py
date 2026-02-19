from Bio import SeqIO
import glob

for file_name in glob.glob("./*.ab1"):
    record = SeqIO.parse(file_name, "abi")
    new_name = file_name.split(".ab1")[0] + "_ab1.fastq"
    count = SeqIO.write(record, new_name, "fastq")
