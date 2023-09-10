import os
import gzip
from sys import argv
myecDNA=argv[1]
tss = "/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/04.ecDNA_detection/ref/tss.bed"
#tss = "/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/04.ecDNA_detection/ref/hg19_tss.bed"
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
with gzip.open(gzfile,"r") as fh:
	line = fh.readline().strip().decode()
	while(line):
		if line.split("\t")[0] in myhashtss:
			s=int(line.split("\t")[1])
			e=int(line.split("\t")[2])
			tss_list=myhashtss[line.split("\t")[0]]
			s_judge_roundregion=[abs(s-x) < 2000 for x in tss_list]
			e_judge_roundregion=[abs(e-x) < 2000 for x in tss_list]
			if True in s_judge_roundregion:
				gene_tss=tss_list[s_judge_roundregion.index(True)]
				ecregion=myhash2[line.split("\t")[0]+":"+str(gene_tss)]
				if s < (gene_tss-24) and e > (gene_tss-24) :
					myhashcell.setdefault(line.split("\t")[3]+":"+ecregion,[]).append("a")
				elif s > (gene_tss-24) and s <(gene_tss+25):
					myhashcell.setdefault(line.split("\t")[3]+":"+ecregion,[]).append("a")
				else:
					myhashcell.setdefault(line.split("\t")[3]+":"+ecregion,[]).append("b")
			elif True in e_judge_roundregion:
				gene_tss=tss_list[e_judge_roundregion.index(True)]
				ecregion=myhash2[line.split("\t")[0]+":"+str(gene_tss)]
				if s < (gene_tss-24) and e > (gene_tss-24) :
					myhashcell.setdefault(line.split("\t")[3]+":"+ecregion,[]).append("a")
				elif s > (gene_tss-24) and s <(gene_tss+25):
					myhashcell.setdefault(line.split("\t")[3]+":"+ecregion,[]).append("a")
				else:
					myhashcell.setdefault(line.split("\t")[3]+":"+ecregion,[]).append("b")
			line = fh.readline().strip().decode()
fh.close()

for i in myhashcell:
	summ= myhashcell[i].count("a")+myhashcell[i].count("b")
	#if myhashcell[i].count("b") != 0:
	num=round(myhashcell[i].count("a")/summ,3)
	out.write(i.split(":")[1]+":"+i.split(":")[2]+"\t"+i.split(":")[0]+"\t"+str(myhashcell[i].count("a"))+"\t"+str(summ)+"\t"+str(num)+"\n")


