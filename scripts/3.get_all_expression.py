#!/usr/bin/python3 
import re
import math
import time
import os
import subprocess
from Bio import SeqIO
from optparse import OptionParser

#################################
#start
 
usage="""
	4.get_all_expresssion.py: estimation of expression level
	Usage: python3 4.get_all_expresssion.py  [-h] [-f miRNA_mature_plant_result.fa] -i in_dir  -o out_dir
"""

parser = OptionParser(usage=usage)
parser.add_option('-i','--in',dest='in_dir',action='store',help='(Required.) '
		+'The path of sRNA files with fasta format.')
parser.add_option('-f','--miRNA_file',dest='miRNA_file',action='store',default="miRNA_database_dir/miRNA_mature_plant_result.fa",help='(Optional.) '
		+'Default: Mature plant miRNA sequences from miRBase for miRNA identification.')
parser.add_option('-o','--out_dir',dest='out_dir',action='store',help='(Required.) '
		+'Output directory of the results.')
(options,args) = parser.parse_args()

if options.in_dir is not None:
    in_path=options.in_dir.rstrip('/')
    filenames = os.listdir(in_path)
else:
    print(parser.print_help())
    exit("Error: The directory of sRNA file is required!")
if options.out_dir is not None:
    out_dir=options.out_dir
    out_dir=out_dir.rstrip('/')
else:
    print(parser.print_help())
    exit("Error: Output directory is required!")

miSeq = options.miRNA_file
start_time = time.time()
print("Run start:")
###################################
def make_dir(filepath):
    if not os.path.exists(filepath):
        os.makedirs(filepath)

###################################
make_dir(out_dir)

mi_dict=dict()
le_dict=dict()
seq_dict=dict()
alcount_dict=dict()
alrpm_dict=dict()
count_dict=dict()
rpm_dict=dict()
title=""

for seq_record in SeqIO.parse(miSeq,"fasta"):
    mi_dict[str(seq_record.seq)] = str(seq_record.id)
    le_dict[str(seq_record.id)] = len(seq_record.seq)

for filename in filenames:
    input_file = in_path + "/" + filename
    for seq_record in SeqIO.parse(input_file,"fasta"):
        if str(seq_record.seq) not in seq_dict:
            seq_dict[seq_record.seq]=1
        else:
            continue

for filename in filenames:
    input_file = in_path + "/" + filename
    title = title + "\t" + filename
    for seq_record in SeqIO.parse(input_file,"fasta"):
        srna_id = str(seq_record.id)
        count_array = srna_id.split("_x")
        rpm_array = srna_id.split("_")
        count_dict[str(seq_record.seq)] = count_array[1]
        rpm_dict[str(seq_record.seq)] = rpm_array[2]

    for key in seq_dict.keys():
        if key in count_dict.keys():
            alcount_dict.setdefault(key,[]).append(count_dict[key])
            alrpm_dict.setdefault(key,[]).append(rpm_dict[key])
        else:
            value2 = 0
            alcount_dict.setdefault(key,[]).append(str(value2))
            alrpm_dict.setdefault(key,[]).append(str(value2))

    count_dict = dict()
    rpm_dict = dict()

outfile = out_dir.split("/")[-1].split("_")
mirna_c_file = out_dir + "/" + outfile[0] + ".miRNA.count.txt"
sirna_c_file = out_dir + "/" + outfile[0] + ".siRNA.count.txt"
srna_c_file = out_dir + "/" + outfile[0] + ".sRNA.count.txt"

mirna_r_file = out_dir + "/" + outfile[0] + ".miRNA.rpm.txt"
sirna_r_file = out_dir + "/" + outfile[0] + ".siRNA.rpm.txt"
srna_r_file = out_dir + "/" + outfile[0] + ".sRNA.rpm.txt"

fm_c = open(mirna_c_file,"w")
fsi_c = open(sirna_c_file,"w")
fs_c = open(srna_c_file,"w")

fm_r = open(mirna_r_file,"w")
fsi_r = open(sirna_r_file,"w")
fs_r = open(srna_r_file,"w")


fm_c.write("sequence" + title + "\n")
fsi_c.write("sequence" + title + "\n")
fs_c.write("sequence" + title + "\n")

fm_r.write("sequence" + title + "\n")
fsi_r.write("sequence" + title + "\n")
fs_r.write("sequence" + title + "\n")

for key in alcount_dict.keys():
    fs_c.write(str(key) + "\t" + "\t".join(alcount_dict[key]) + "\n")
    fs_r.write(str(key) + "\t" + "\t".join(alrpm_dict[key]) + "\n")
    if str(key) in mi_dict.keys():
        fm_c.write(str(key) + "\t" + "\t".join(alcount_dict[key]) + "\n")
        fm_r.write(str(key) + "\t" + "\t".join(alrpm_dict[key]) + "\n")
    else:
        if (len(key)>=21) and (len(key)<=24):
            fsi_c.write(str(key) + "\t" + "\t".join(alcount_dict[key]) + "\n")
            fsi_r.write(str(key) + "\t" + "\t".join(alrpm_dict[key]) + "\n")
        else:
            continue

run_time=int(time.time() - start_time)
print("Run complete [3.get_all_expression.py]: "+"%d seconds elapsed " % run_time)
