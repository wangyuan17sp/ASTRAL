#!/usr/bin/env python
import sys
import os
from optparse import OptionParser
import numpy as np
import itertools
import random
import subprocess
import re
if __name__ == '__main__':
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option("-i","--input",dest="inpt",type="string",
			help="read pool of branches from FILENAME")
	parser.add_option("-s","--species",dest="sp",type="string",
			help="read species info from FILENAME")
	parser.add_option("-o","--output",dest="out",type="string",
			help="the PATH to write the generated files")
	parser.add_option("-v","--verbose",dest="verbose",
			help="Verbose",default=1)
	(options,args) = parser.parse_args()
	inpt = options.inpt
	outpath = options.out
	sp = options.sp
	if ( not options.inpt  or not options.out or not options.sp):
		sys.exit("Please enter pool of brnaches file, species info file, and output folder location")

	poolSpeciesBranches = dict()

	branches = dict()
	f = open(sp, 'r')
	for x in f:
		y = re.search('^\{[0-9]',x)
		if y:
			bpInfo=re.sub('.*\[','[',x)
			bpInfo=re.sub('[\[\]{} ]','',bpInfo)
			bpInfo=re.sub(':.*','',bpInfo)
			poolSpeciesBranches[bpInfo] = 1;
	f.close()
	g = open(inpt,'r')
	for line in g:
		y = re.search(':',line)
		if y:
			bpInfo = re.sub(':.*','',line)
			pp     = re.sub('.*:','',line)
			if bpInfo in poolSpeciesBranches:
				branches[bpInfo] = str(1)+" "+str(pp)
			else:
				branches[bpInfo] = str(0)+" "+str(pp)
	g.close()
	f = open(outpath+"/ppOfBranches.txt",'w')
	f.close()
	f = open(outpath+"/ppOfBranches.txt",'w')
	for key in branches:
		print >> f,branches[key].replace("\n","")
	f.close()
