import os
from sys import argv
import math

candidate_bin=argv[1]
candidate_pos=argv[2]
outfile=argv[3]
myhash={}
out=open(outfile,"wt")

with open (candidate_bin,"r") as fh:
	line=fh.readline().strip()
	while(line):
		myhash.setdefault(line.split("\t")[0],[]).append(line.split("\t")[1]) 
		line=fh.readline().strip()
fh.close()

with open (candidate_pos,"r") as fh:
	line= fh.readline().strip()
	chrr =line.split("\t")[0]
	start = line.split("\t")[1]
	end = line.split("\t")[2]
	line= fh.readline().strip()
	tmp = chrr + ":" + start + "_" + end
	while(line):
		newchrr =line.split("\t")[0]
		newstart = line.split("\t")[1]
		newend = line.split("\t")[2]
		tmp1 = newchrr + ":" + newstart + "_" + newend
		if newchrr==chrr and newstart==end:
			if len(list(set(myhash[tmp]) & set(myhash[tmp1]))) > 7 :
				end = newend
			else :
				out.write(chrr+"\t"+start+"\t"+end+"\n")
				chrr=newchrr
				start=newstart
				end=newend
		else :
			out.write(chrr+"\t"+start+"\t"+end+"\n")
			chrr=newchrr
			start=newstart
			end=newend
		tmp = tmp1
		line= fh.readline().strip()
