import os
from sys import argv
import re
junctionsite=argv[1]
region=argv[2]
outfile=argv[3]
out=open(outfile,"wt")
myhash={}

with open (junctionsite,"r") as fh1:
	line=fh1.readline().strip()
	while(line):
		chrom = line.split("\t")[1]
		sseq1_s = line.split("\t")[2]
		sseq1_e = line.split("\t")[3]
		sseq2_s = line.split("\t")[6]
		sseq2_e = line.split("\t")[7]
		myhash.setdefault(chrom,[]).append(sseq1_s+":"+sseq1_e+"_"+sseq2_s+":"+sseq2_e+"_"+line.split("\t")[13])
		line=fh1.readline().strip()

with open (region,"r") as fh2:
	fh2.readline().strip()
	line=fh2.readline().strip()
	while(line):
		chrr=line.split("\t")[0].split(":")[0]
		my_s=int(line.split("\t")[0].split(":")[1].split("_")[0])
		my_e=int(line.split("\t")[0].split(":")[1].split("_")[1])
		if chrr in myhash :
			for i in myhash[chrr]:
				sseq1_s=int(i.split("_")[0].split(":")[0])
				sseq1_e=int(i.split("_")[0].split(":")[1])
				sseq2_s=int(i.split("_")[1].split(":")[0])
				sseq2_e=int(i.split("_")[1].split(":")[1])
				if sseq1_s > my_s and sseq1_s < my_e or sseq1_e > my_s and sseq1_e < my_e :
					out.write("oneseq"+"\t"+line+"\t"+i.split("_")[0]+"\t"+i.split("_")[1]+"\t"+i.split("_")[2]+"\n")
					if sseq2_s > my_s and sseq2_s < my_e or sseq2_e > my_s and sseq2_e < my_e :
						out.write("twoseq"+"\t"+line+"\t"+i.split("_")[0]+"\t"+i.split("_")[1]+"\t"+i.split("_")[2]+"\n")
				elif sseq2_s > my_s and sseq2_s < my_e or sseq2_e > my_s and sseq2_e < my_e :
					out.write("oneseq"+"\t"+line+"\t"+i.split("_")[0]+"\t"+i.split("_")[1]+"\t"+i.split("_")[2]+"\n")
				elif sseq1_s > my_s-500000 and sseq1_s < my_e+500000 or sseq1_e > my_s-500000 and sseq1_e < my_e+500000 :
					out.write("neighborhood"+"\t"+line+"\t"+i.split("_")[0]+"\t"+i.split("_")[1]+"\t"+i.split("_")[2]+"\n")
				elif sseq2_s > my_s-500000 and sseq2_s < my_e+500000 or sseq2_e > my_s-500000 and sseq2_e < my_e+500000 :
					out.write("neighborhood"+"\t"+line+"\t"+i.split("_")[0]+"\t"+i.split("_")[1]+"\t"+i.split("_")[2]+"\n")
		line=fh2.readline().strip()		
