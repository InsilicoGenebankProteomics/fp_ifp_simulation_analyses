import matplotlib.pyplot as plt
import numpy as np
import sys

nfiles=len(sys.argv)-1

micradius=[]
for j in range(1,int(nfiles/2)+1):
	cmdistancesfile=str(sys.argv[j])
	print(cmdistancesfile)
	file = open(cmdistancesfile)
	lines = file.readlines()
	file.close()

	nlines = len(lines)
	nvars = len(lines[2].strip().split())
	print(nlines,nvars)
	for i in range(2,nlines):
		l=lines[i].strip().split()
		micradius.append(float(l[1]))
plt.hist(micradius,bins=100,range=(22,32),alpha=0.5,label='dpc')

print('Second batch')
micradius=[]
for j in range(int(nfiles/2)+1,nfiles+1):
	cmdistancesfile=str(sys.argv[j])
	print(cmdistancesfile)
	file = open(cmdistancesfile)
	lines = file.readlines()
	file.close()

	nlines = len(lines)
	nvars = len(lines[2].strip().split())
	print(nlines,nvars)
	for i in range(2,nlines):
		l=lines[i].strip().split()
		micradius.append(float(l[1]))
plt.hist(micradius,bins=100,alpha=0.5,range=(22,32),label='lmpg')

plt.legend()
plt.show()