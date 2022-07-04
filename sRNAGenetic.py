#/usr/bin/python3
import re
import math
import time
import os
import subprocess
from Bio import SeqIO
from optparse import OptionParser

usage="""
    sRNAGenetic.py: miRNA identification and calculate count and rpm of sRNA
    Usage: python3 sRNAGenetic.py [-h] [--ft] [--mt 15] [--lt 30] [-f miRNA_mature_plant_result.fa] -i in_dir  -o out_dir
"""
parser = OptionParser(usage=usage)
parser.add_option('-i','--in',dest='in_dir',action='store',help='(Required.) '
                 +'The path of sRNA files with fasta format.')
parser.add_option('--ft','--file_type',dest='file_type',action='store',help='(Required.) '
                 +'The format of input file. The option "SRA" represents the clean sRNA fasta-format file after trimming and QC.'
                 +'The option "GEO" represents the txt files of sRNA data that download from GEO database.')
parser.add_option('--mt','--more_than',dest='more_than',action='store',default=15,help='(Optional.) '
                 +'Threoshold of sRNA minimum length. sRNA with length less than MORE_THAN will be '
                 +'classified as non-sRNA. The Default is 15 .')
parser.add_option('--lt','--less_than',dest='less_than',action='store',default=30,help='(Optional.) '
                 +'Threoshold of sRNA maximum length. sRNA with length more than LESS_THAN will be '
                 +'classified as non-sRNA. The Default is 30 .')
parser.add_option('-f','--miRNA_file',dest='miRNA_file',action='store',default="miRNA_database_dir/miRNA_mature_plant_result.fa",help='(Optional.) '
                +'Default: miRNA_database_dir/miRNA_mature_plant_result.fa: Mature plant miRNA sequences from miRBase for miRNA identification.')
parser.add_option('-o','--out_dir',dest='out_dir',action='store',help='(Required.) '
                 +'Output directory of the results.')
(options,args) = parser.parse_args()

if (options.file_type == "SRA") or (options.file_type == "GEO"):
    FT = options.file_type
else:
    exit("Error: Please select the input file type correctly!") 
if options.in_dir is not None:
    in_path=options.in_dir.rstrip('/')
    dirs = os.listdir(in_path)
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
miSeq = options.miRNA_file
start_time = time.time()
print("Run start:")

#################################
def make_dir(filepath):
    if not os.path.exists(filepath):
        os.makedirs(filepath)

#################################
make_dir(out_dir)

#step1:Triming && QC
#step2:sRNA filtration
for file_dir in dirs:
	print("For: "+file_dir+"\n"+"step1:sRNA filtration")
	step1_in_path = in_path + "/" + file_dir + "/"
	step1_out_dir = out_dir + "/" + file_dir + "/sRNA_dir/"
	if (str(FT) == "SRA"):
		subprocess.run(["python3","scripts/1.sRNA_fasta_filtration.py","--mt",str(MT_len),"--lt",str(LT_len),"-i",step1_in_path,"-o",step1_out_dir])
	elif (str(FT) == "GEO"):
		subprocess.run(["python3","scripts/2.sRNA_txt_filtration.py","--mt",str(MT_len),"--lt",str(LT_len),"-i",step1_in_path,"-o",step1_out_dir])
	else:
		exit("Error: Please select the input file type correctly!")

#step3:expression level
for file_dir in dirs:
	print("For: "+file_dir+"\n"+"step2:Caculation of expression level (miRNA,siRNA,sRNA)")
	step2_in_dir = out_dir + "/" + file_dir + "/sRNA_dir/"
	step2_out_dir = out_dir + "/" + file_dir + "/"
	subprocess.run(["python3","scripts/3.get_all_expression.py","-f",miSeq,"-i",step2_in_dir,"-o",step2_out_dir])

run_time=int(time.time() - start_time)
min_time=round(run_time/60,0)
print("Run complete: "+"%d minutes elapsed " % min_time)

