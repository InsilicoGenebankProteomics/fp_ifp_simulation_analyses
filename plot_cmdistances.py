import matplotlib.pyplot as plt
import numpy as np
import sys

cmdistancesfile=str(sys.argv[1])
file = open(cmdistancesfile)
lines = file.readlines()
file.close()

nlines = len(lines)
nvars = len(lines[2].strip().split())
rescm = np.zeros((nlines-2,nvars-4))
print(nlines,nvars)
time=[]
micradius=[]
micradius_var=[]
pepcm=[]
for i in range(2,nlines):
	l=lines[i].strip().split()
	time.append((i-2)*0.01)
	micradius.append(float(l[1]))
	micradius_var.append(float(l[2]))
	pepcm.append(float(l[3]))
	rescm[i-2,:]=l[4:]

micradius=np.array(micradius)
micradius_var=np.array(micradius_var)
print(micradius,micradius_var)	
plt.plot(time,micradius)#,label='micelle radius')
plt.fill_between(time, micradius-micradius_var, micradius+micradius_var, alpha=0.5)
#plt.plot(time,pepcm)
plt.plot(time,rescm[:,0],label='Cterminal')
plt.plot(time,rescm[:,-1],label='Nterminal')
plt.plot(time,np.mean(rescm,axis=1),label='peptide center')
plt.ylim(0,35)
plt.xlim(0,150)
plt.legend()
plt.show()