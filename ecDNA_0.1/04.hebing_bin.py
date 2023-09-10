import os
from sys import argv
inputfile = argv[1]
outfile=argv[2]
myhash={}
out=open(outfile,"wt")
with open (inputfile,"r") as fh:
	line= fh.readline().strip()
	chrr =line.split("\t")[0].split(":")[0]
	start = line.split("\t")[0].split(":")[1].split("_")[0]
	end = line.split("\t")[0].split(":")[1].split("_")[1]
	line= fh.readline().strip()
	while(line):
		newchrr =line.split("\t")[0].split(":")[0]
		newstart = line.split("\t")[0].split(":")[1].split("_")[0]
		newend = line.split("\t")[0].split(":")[1].split("_")[1]
		if newchrr==chrr and newstart==end:
			end = newend
		else :
			out.write(chrr+"\t"+start+"\t"+end+"\n")
			chrr=newchrr
			start=newstart
			end=newend
		line= fh.readline().strip()

