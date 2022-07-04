#/usr/bin/python3
import re
import math
import time
import os
from Bio import SeqIO
from optparse import OptionParser

usage="""
    2.sRNA_txt_filtration.py: select the Reads with appropriate length and calculate the RPM and Count of sRNA
    Usage: python3  2.sRNA_txt_filtration.py [-h] [--mt 15] [--lt 30] -i in_dir  -o out_dir
"""
parser = OptionParser(usage=usage)
parser.add_option('-i','--in',dest='in_dir',action='store',help='(Required.) '
                 +'The path of sRNA files with fasta format.')
parser.add_option('--mt','--more_than',dest='more_than',action='store',default=15,help='(Optional.) '
		 +'Threoshold of sRNA minimum length. sRNA with length less than MORE_THAN will be '
                 +'classified as non-sRNA. The Default is 15 .')
parser.add_option('--lt','--less_than',dest='less_than',action='store',default=30,help='(Optional.) '
                 +'Threoshold of sRNA maximum length. sRNA with length more than LESS_THAN will be '
                 +'classified as non-sRNA. The Default is 30 .')
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

MT_len = int(options.more_than)
LT_len = int(options.less_than)
start_time = time.time()
print("Run start:")
#################################
def make_dir(filepath):
    if not os.path.exists(filepath):
        os.makedirs(filepath)

#################################
make_dir(out_dir)
seq_dict = dict ()
L_num = 0
total = 0
number = 0
All_num = 0
Ratio = 0

for filename in filenames:
	print("For sRNA file:" + str(filename))
	input_file = in_path + "/" + filename
	with open(input_file) as txt_object:
		for txt_line in txt_object:
			txt_line = txt_line.rstrip("\n")
			txt_array = txt_line.split("\t")
			if (re.findall(r"[^ATCGatcg]",str(txt_array[0]))):
				continue
			else:
				All_num = All_num + int(txt_array[1])
				if (len(txt_array[0])>=MT_len) and (len(txt_array[0])<=LT_len):
					seq_dict[txt_array[0]] = int(txt_array[1])
					total = total + int(txt_array[1])
				else:
					L_num = L_num + int(txt_array[1])
	outfile = filename.split(".")[0]
	output_file = out_dir + "/" + outfile + ".filiter.fasta"
	fd = open(output_file,"w")
	for key in seq_dict.keys():
		number += 1
		rpm = round(seq_dict[key]/total*1000000,2)
		fd.write(">"+str(outfile)+"_"+str(number)+"_"+str(rpm)+"_x"+str(seq_dict[key])+"\n")
		fd.write(str(key)+"\n")

	Ratio = round((total/All_num)*100,2)
	print("Total Number of Reads:" + str(All_num))
	print("Number of reads with length limitation filtered out:" + str(L_num))
	#print("Number of reads with N base filtered out:" + str(N_num))
	print("Number of valid reads(ratio):" + str(total) + "(" + str(Ratio),"%)")
	seq_dict = dict ()
	L_num = 0
	total = 0
	number = 0
	All_num = 0
	Ratio = 0

run_time=int(time.time() - start_time)
print("Run complete [2.sRNA_txt_filtration.py]: "+"%d seconds elapsed " % run_time)
