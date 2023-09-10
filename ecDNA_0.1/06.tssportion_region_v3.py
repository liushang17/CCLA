import os
import gzip
from sys import argv
import numpy as np
myecDNA=argv[1]
#tss = "/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/04.ecDNA_detection/ref/tss.bed"
#tss = "/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/04.ecDNA_detection/ref/hg19_tss.bed"
tss = "/hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/ref/mm10/tss.bed"
gzfile=argv[2]
outfile=argv[3]
myhash={}
out=open(outfile,"wt")

with open (myecDNA,"r") as fh:
	line= fh.readline().strip()
	while (line):
		myhash[line.split("\t")[1]]=line.split("\t")[0]
		line= fh.readline().strip()
fh.close()
myhashtss={}
myhash2={}
with open (tss,"r") as fh:
	line = fh.readline().strip()
	while (line):
		if line.split("\t")[3] in myhash:
			region=myhash[line.split("\t")[3]]
			chrr=line.split("\t")[0]
			tss=line.split("\t")[1]
			myhashtss.setdefault(chrr,[]).append(int(tss))
			myhash2[chrr+":"+tss]=region
		line= fh.readline().strip()

myhashcell={}
with open(gzfile,"r") as fh:
	line = fh.readline().strip()
	while(line):
		if line.split("\t")[0] in myhashtss:
			s=int(line.split("\t")[1])
			e=int(line.split("\t")[2])
			tss_list=np.array(myhashtss[line.split("\t")[0]])
			tss_list1=tss_list-24
			tss_list2=tss_list+25			
			s_judge_roundregion=(np.where(abs(tss_list - s) < 2000)[0])
			e_judge_roundregion=(np.where(abs(tss_list - e) < 2000)[0])
			judge_roundregion = len(s_judge_roundregion) + len(e_judge_roundregion)

			judge_roundtss1=len(np.where(tss_list1-e > 0)[0])
			judge_roundtss2=len(np.where(tss_list2-s < 0)[0])
			judge_roundtss=len(tss_list) - (judge_roundtss1 + judge_roundtss2)

			if judge_roundregion > 0:
				gene_tss1=tss_list[s_judge_roundregion].tolist()
				gene_tss2=tss_list[e_judge_roundregion].tolist()
				gene_tss = gene_tss1 + gene_tss2
				ecregion=myhash2[line.split("\t")[0]+":"+str(gene_tss[0])]
				if judge_roundtss > 0:
					myhashcell.setdefault(line.split("\t")[3]+":"+ecregion,[]).append("a")
				else:
					myhashcell.setdefault(line.split("\t")[3]+":"+ecregion,[]).append("b")
		line = fh.readline().strip()
fh.close()

for i in myhashcell:
	summ= myhashcell[i].count("a")+myhashcell[i].count("b")
	#if myhashcell[i].count("b") != 0:
	num=round(myhashcell[i].count("a")/summ,3)
	out.write(i.split(":")[1]+":"+i.split(":")[2]+"\t"+i.split(":")[0]+"\t"+str(myhashcell[i].count("a"))+"\t"+str(summ)+"\t"+str(num)+"\n")


