import glob

with open("all_mosaic_clones.fastq", "w") as output:
    for file_name in glob.glob("C*.fastq"):
        with open(file_name, "r") as entry:
            i = 0
            for line in entry:
                if i == 0:
                    output.write("@" + file_name.split(".fastq")[0] + "\n")
                    i += 1
                else:
                    output.write(line.strip() + "\n")
                print(line)
