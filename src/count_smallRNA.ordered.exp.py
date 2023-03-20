import sys 
sys.path.insert(0, './')
#import dr_tools, pysam
import argparse, os


parser = argparse.ArgumentParser()
parser.add_argument('-a', '--anno', help='reads order',required=True)
parser.add_argument('-c', '--tmpOvcount', required=True, help='counts overlap different mature ncRNA')
parser.add_argument('-o', '--output', required=True, help='output directory')
o = parser.parse_args()


type_order={}
with open(o.anno,mode="r") as file:
    for line in file:
        line=line.strip("\n").split("\t")
        type_order[line[0]]=int(line[1])


reads_brief={}
reads_type={}
with open(o.tmpOvcount,mode="r") as file:
    for line in file:
        line=line.strip("\n").split("\t")
        if line[2] not in type_order:
#            continue
            line[2]="OTHER"
        if line[0] in reads_brief:
            if type_order[line[2]] > type_order[reads_type[line[0]]]:
                continue
            elif type_order[line[2]] == type_order[reads_type[line[0]]]:
                if line[1] not in reads_brief[line[0]]:
                    reads_brief[line[0]].append(line[1])
            elif type_order[line[2]] < type_order[reads_type[line[0]]]:
                reads_brief[line[0]]=[]
                reads_brief[line[0]].append(line[1])
                reads_type[line[0]]=line[2]

                    
        else:
            reads_brief[line[0]]=[]
            reads_brief[line[0]].append(line[1])
            reads_type[line[0]]=line[2]
with open(o.output+".od.exp.tmp.txt","w") as output_exp_tmp:
    anno_count={}
    for read in reads_brief:
        for anno in reads_brief[read]:
            if anno not in anno_count:
                anno_count[anno]=1/len(reads_brief[read])
            else:
                anno_count[anno]=1/len(reads_brief[read])+ anno_count[anno]
            print(read,anno,reads_type[read],sep="\t",file=output_exp_tmp)

output_exp_tmp.close()

with open(o.output+".od.exp.txt","w") as output_exp:
    for anno in  anno_count:
        print(anno,anno_count[anno],sep="\t",file=output_exp)
        
output_exp.close()       
