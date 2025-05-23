import matplotlib.pyplot as plt
import numpy as np
import sys

nfiles=len(sys.argv)-1
for j in range(nfiles):
	interactionsfile=str(sys.argv[1])
	file = open(interactionsfile)
	lines = file.readlines()
	file.close()

	perres=True
	res = []
	water = []
	lipid = []
	i=1
	l=lines[i].strip().split()
	while perres:
		res.append(l[1])
		water.append(float(l[2]))
		lipid.append(float(l[3]))
		print(l[1:4])
		i=i+1
		l=lines[i].strip().split()
		if l[0] == '#Interactions': perres=False
	if j==0:
		water_ave=np.array(water)
		lipid_ave=np.array(lipid)
	else:
		water_ave=water_ave+np.array(water)
		lipid_ave=lipid_ave+np.array(lipid)

water_ave=water_ave/nfiles
lipid_ave=lipid_ave/nfiles
plt.bar(res, water, label='water')
plt.bar(res, lipid, bottom=water, label='lipid')
plt.legend()
plt.show()