#!/usr/bin/python3
from optparse import OptionParser
import os
#################################
#start
 
usage="""
	Usage: python3 get_merged_filiter_by_count.py  --replicate 3 --count 5 --P1 P1.siRNA.count.txt  --P2 P2.siRNA.count.txt  --F1 F1.siRNA.count.txt -o out_dir/
"""

parser = OptionParser(usage=usage)
parser.add_option('--replicate','--replicate',dest='replicate',action='store',help='(Required.) '
                +'Number of biological replicates.')
parser.add_option('--count','--count',dest='count',action='store',default=5,help='(Optional.) '
                +'Default: Filiter criteria for count value.')
parser.add_option('--P1','--P1',dest='P1',action='store',help='(Required.) '
		+'Count table of P1 species.')
parser.add_option('--P2','--P2',dest='P2',action='store',help='(Required.) '
		+'Count table of P2 species.')
parser.add_option('--F1','--F1',dest='F1',action='store',help='(Required.) '
		+'Count table of F1 species.')
parser.add_option('-o','--out_dir',dest='out_dir',action='store',help='(Required.) '
                +'Output directory of the results.')
(options,args) = parser.parse_args()

cova = int(options.count)

if options.P1 is not None:
    P1_id =  options.P1
else:
    print(parser.print_help())
    exit("Error: The count file of 'P1' is required!")

if options.P2 is not None:
    P2_id =  options.P2
else:
    print(parser.print_help())
    exit("Error: The count file of 'P2' is required!")

if options.F1 is not None:
    F1_id =  options.F1
else:
    print(parser.print_help())
    exit("Error: The count file of 'F1' is required!")

if options.replicate is not None:
    re_nu =  options.replicate
else:
    print(parser.print_help())
    exit("Error: The number of biological replicates is required!")

if options.out_dir is not None:
    out_dir=options.out_dir
    out_dir=out_dir.rstrip('/')
else:
    print(parser.print_help())
    exit("Error: Output directory is required!")

###################################
def make_dir(filepath):
    if not os.path.exists(filepath):
        os.makedirs(filepath)

def count_filiter(list1,cova):
    for n in list1:
        n = int(n)
        if n < cova:
            return("False")
        else:
            continue
    return("True")
###################################
make_dir(out_dir)

###################################
i = 0
A_dict = dict()
B_dict = dict()
C_dict = dict()
id_dict = dict()

x = ["0"]*int(re_nu)
y = "\t".join(x)

with open(P1_id) as file1:
    for Aline in file1:
        id1,count1 = Aline.rstrip("\n").split("\t",1)
        A_dict[id1] = str(count1)
        id_dict[id1] = 1

with open(P2_id) as file2:
    for Bline in file2:
        id2,count2 = Bline.rstrip("\n").split("\t",1)
        B_dict[id2] = str(count2)
        id_dict[id2] = 1

with open(F1_id) as file3:
    for Cline in file3:
        id3,count3 = Cline.rstrip("\n").split("\t",1)
        C_dict[id3] = str(count3)
        id_dict[id3] = 1

a = 0
b = 0
c = 0
for key in id_dict.keys():
    if key not in A_dict.keys():
        A_dict[key] = y
    else:
        continue
        #a = a + 1

for key in id_dict.keys():
    if key not in B_dict.keys():
        B_dict[key] = y
    else:
        continue
        #b = b + 1

for key in id_dict.keys():
    if key not in C_dict.keys():
        C_dict[key] = y
    else:
        continue
        #c = c + 1

P1_count_file = out_dir + "/" + "P1.siRNA.filiter.count.txt"
P2_count_file = out_dir + "/" + "P2.siRNA.filiter.count.txt"
F1_count_file = out_dir + "/" + "F1.siRNA.filiter.count.txt"
merge_count = out_dir + "/" + "merged.filiter.count.txt"

P1_c = open(P1_count_file,"w")
P2_c = open(P2_count_file,"w")
F1_c = open(F1_count_file,"w")
me_c = open(merge_count,"w")

for key in id_dict.keys():
    A_list = A_dict[key].split("\t")
    B_list = B_dict[key].split("\t")
    C_list = C_dict[key].split("\t")
    Alo = count_filiter(A_list,cova)
    Blo = count_filiter(B_list,cova)
    Clo = count_filiter(C_list,cova)
    if Alo == "True" or Blo == "True" or Clo == "True":
        P1_c.write(key + "\t" + A_dict[key] + "\n")
        P2_c.write(key + "\t" + B_dict[key] + "\n")
        F1_c.write(key + "\t" + C_dict[key] + "\n")
        me_c.write(key + "\t" + A_dict[key] + "\t" + B_dict[key] + "\t" + C_dict[key] + "\n")
    else:
        continue

