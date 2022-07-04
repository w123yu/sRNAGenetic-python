#!/usr/bin/python3
from optparse import OptionParser

#################################
#start
 
usage="""
	Usage: python3 orthogroups.py  --F1P1_vs_P1 F1P1.vs.P1.id.txt  --F1P2_vs_P2 F1P2.vs.P2.id.txt  --F1P1_vs_F1P2 F1P1.vs.F1P2.id.txt > results.txt
"""

parser = OptionParser(usage=usage)
parser.add_option('--F1P1_vs_P1','--F1P1_vs_P1',dest='F1P1_vs_P1',action='store',help='(Required.) '
		+'Reciprocal best gene pairs between the genome of P1 species and one sub-genome in F1 that are of the same type as the parental P1 genome.')
parser.add_option('--F1P2_vs_P2','--F1P2_vs_P2',dest='F1P2_vs_P2',action='store',help='(Required.) '
		+'Reciprocal best gene pairs between the genome of P2 species and one sub-genome in F1 that are of the same type as the parental P2 genome.')
parser.add_option('--F1P1_vs_F1P2','--F1P1_vs_F1P2',dest='F1P1_vs_F1P2',action='store',help='(Required.) '
		+'Reciprocal best gene pairs between sub-genomes of F1.')
(options,args) = parser.parse_args()

if options.F1P1_vs_P1 is not None:
    F1P1_vs_P1_id =  options.F1P1_vs_P1
else:
    print(parser.print_help())
    exit("Error: The file of 'F1P1_vs_P1' is required!")

if options.F1P2_vs_P2 is not None:
    F1P2_vs_P2_id =  options.F1P2_vs_P2
else:
    print(parser.print_help())
    exit("Error: The file of 'F1P2_vs_P2' is required!")

if options.F1P1_vs_F1P2 is not None:
    F1P1_vs_F1P2_id =  options.F1P1_vs_F1P2
else:
    print(parser.print_help())
    exit("Error: The file of 'F1P1_vs_F1P2' is required!")

###################################
i = 0
A_dict = dict()
B_dict = dict()
C_dict = dict()

with open(F1P1_vs_P1_id) as file1:
    for Aline in file1:
        array1 = Aline.rstrip("\n").split("\t")
        A_dict[str(array1[0])] = str(array1[1])

with open(F1P2_vs_P2_id) as file2:
    for Bline in file2:
        array2 = Bline.rstrip("\n").split("\t")
        B_dict[str(array2[1])] = str(array2[0])

with open(F1P1_vs_F1P2_id) as file3:
    for Cline in file3:
        array3 = Cline.rstrip("\n").split("\t")
        C_dict[str(array3[1])] = str(array3[0])

for key in A_dict.keys():
    Bnapu_At = A_dict[key]
    if Bnapu_At in C_dict.keys():
        Bnapu_Ct = C_dict[Bnapu_At]
        if Bnapu_Ct in B_dict.keys():
            i = i + 1
            print("homoeolog" + str(i) + "\t" + key + "\t" + str(B_dict[Bnapu_Ct]) + "\t" + Bnapu_At + "\t" + Bnapu_Ct)


