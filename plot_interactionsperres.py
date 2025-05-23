import matplotlib.pyplot as plt
import numpy as np
import sys

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
	
plt.bar(res, water)
plt.bar(res, lipid, bottom=water)
plt.show()