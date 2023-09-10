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
with open (tss,"r") as fh:
	line = fh.readline().strip()
	while (line):
		if line.split("\t")[3] in myhash:
			chrr=line.split("\t")[0]
			tss=int(line.split("\t")[1])
			tssregion=str(tss-24)+"_"+str(tss+25)
			roundregion=str(tss-24-2000)+"_"+str(tss+25+2000)
			myhashtss.setdefault(chrr,[]).append(line.split("\t")[3]+":"+tssregion+":"+roundregion)
		line= fh.readline().strip()
fh.close()
myhashcell={}
with gzip.open(gzfile,"r") as fh:
	line = fh.readline().strip().decode()
	while(line):
		if line.split("\t")[0] in myhashtss	:
			for i in myhashtss[line.split("\t")[0]]:
				tsss= int(i.split(":")[1].split("_")[0])
				tesse=int(i.split(":")[1].split("_")[1])
				rounds=int(i.split(":")[2].split("_")[0])
				rounde=int(i.split(":")[2].split("_")[1])
				if int(line.split("\t")[1]) > rounds and int(line.split("\t")[1]) < rounde:
						region=myhash[i.split(":")[0]]
						if  int(line.split("\t")[1]) < tesse and  int(line.split("\t")[1]) > tsss:
							myhashcell.setdefault(line.split("\t")[3]+":"+region,[]).append("a")
						elif int(line.split("\t")[2]) < tesse and  int(line.split("\t")[2]) > tsss:
							myhashcell.setdefault(line.split("\t")[3]+":"+region,[]).append("a")
						else:
							myhashcell.setdefault(line.split("\t")[3]+":"+region,[]).append("b")
						break
				elif int(line.split("\t")[2]) > rounds and int(line.split("\t")[2]) <rounde:
						region=myhash[i.split(":")[0]]
						if  int(line.split("\t")[1]) < tesse and  int(line.split("\t")[1]) > tsss:
							myhashcell.setdefault(line.split("\t")[3]+":"+region,[]).append("a")
						elif int(line.split("\t")[2]) < tesse and  int(line.split("\t")[2]) > tsss:
							myhashcell.setdefault(line.split("\t")[3]+":"+region,[]).append("a")
						else:
							myhashcell.setdefault(line.split("\t")[3]+":"+region,[]).append("b")
						break
		line = fh.readline().strip().decode()
fh.close()
for i in myhashcell:
	summ= myhashcell[i].count("a")+myhashcell[i].count("b")
	#if myhashcell[i].count("b") != 0:
	num=round(myhashcell[i].count("a")/summ,3)
	out.write(i.split(":")[1]+":"+i.split(":")[2]+"\t"+i.split(":")[0]+"\t"+str(myhashcell[i].count("a"))+"\t"+str(myhashcell[i].count("b"))+"\t"+str(num)+"\n")
