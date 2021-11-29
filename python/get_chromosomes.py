import os
import sys
import HTSeq

def get_chromosomes(chrs_in_file_path, out_file_path):
    wd = os.getcwd()
    out_file = open(wd + "/" + out_file_path, "w")
    chrs = [(seq.name, seq) for seq in HTSeq.FastaReader(wd + "/" + chrs_in_file_path)]
    for (name, seq) in chrs:
        out_file.write(name + "\t" + str(len(seq)) + "\n")
    out_file.close()



if len(sys.argv) < 3:
    print("Too few arguments. First argument: chrs_in_file_path, Second argument: out_file_path")
    sys.exit()

get_chromosomes(sys.argv[1], sys.argv[2])
