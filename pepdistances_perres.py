import matplotlib.pyplot as plt
import numpy as np
import sys

nfiles=len(sys.argv)-1
meanpepcm=[]
meanrescm=[]
stdrescm=[]
meanmicelleradius=[]
lipstdmicelleradius=[]
for j in range(1,int(nfiles/2)+1):
	cmdistancesfile=str(sys.argv[j])
	print(cmdistancesfile)
	file = open(cmdistancesfile)
	lines = file.readlines()
	file.close()

	nlines = len(lines)
	nvars = len(lines[2].strip().split())
	rescm = np.zeros((nlines-2,nvars-4))
	print(nlines,nvars)
	for i in range(2,nlines):
		l=lines[i].strip().split()
		rescm[i-2,:]=l[4:]
		meanmicelleradius.append(float(l[1]))
		lipstdmicelleradius.append(float(l[2]))
	pepcm=np.mean(rescm,axis=1)
	meanpepcm.append(np.mean(pepcm))
	meanrescm.append(np.mean(rescm,axis=0))
	stdrescm.append(np.std(rescm,axis=0))
print(rescm)
meanmicelleradius=np.mean(meanmicelleradius)
lipstdmicelleradius=np.mean(lipstdmicelleradius)
meanpepcm=np.mean(meanpepcm,axis=0)
meanrescm=np.mean(meanrescm,axis=0)
stdrescm=np.mean(stdrescm,axis=0)
print(meanrescm)
xres=np.arange(1,len(meanrescm)+1)
pepcolor='#06d6a0'
miccolor='#ffd166'

plt.plot(xres,meanrescm,color=pepcolor,label='peptide')
plt.fill_between(xres, meanrescm-stdrescm, meanrescm+stdrescm,color=pepcolor, alpha=0.5)
plt.axhline(y=meanpepcm, linestyle='--',color=pepcolor)
plt.axhline(y=meanmicelleradius, color=miccolor, linestyle='--',label='micelle')
plt.fill_between(xres, meanmicelleradius-lipstdmicelleradius, meanmicelleradius+lipstdmicelleradius,color=miccolor, alpha=0.5)
plt.legend()

plt.title('ifp-dpc')
resnames=['1G','2A','3A','4L','5Q','6I','7P','8F','9A','10M','11Q','12M','13A','14Y','15R','16F']
#resnames=['1M','2W','3K','4T','5P','6T','7L','8K','9Y','10F','11G','12G','13F','14N','15F','16S','17Q','18I','19L']
plt.xticks(xres,resnames,fontsize=8, rotation=60)
plt.ylabel("Distance (Å)")
plt.ylim(0,35)
plt.xlim(1,len(xres))
plt.show()

print('Second batch')
meanpepcm=[]
meanrescm=[]
stdrescm=[]
meanmicelleradius=[]
lipstdmicelleradius=[]
for j in range(int(nfiles/2)+1,nfiles+1):
	cmdistancesfile=str(sys.argv[j])
	print(cmdistancesfile)
	file = open(cmdistancesfile)
	lines = file.readlines()
	file.close()

	nlines = len(lines)
	nvars = len(lines[2].strip().split())
	rescm = np.zeros((nlines-2,nvars-4))
	print(nlines,nvars)
	for i in range(2,nlines):
		l=lines[i].strip().split()
		rescm[i-2,:]=l[4:]
		meanmicelleradius.append(float(l[1]))
		lipstdmicelleradius.append(float(l[2]))
	pepcm=np.mean(rescm,axis=1)
	meanpepcm.append(np.mean(pepcm))
	meanrescm.append(np.mean(rescm,axis=0))
	stdrescm.append(np.std(rescm,axis=0))
print(rescm)
meanmicelleradius=np.mean(meanmicelleradius)
lipstdmicelleradius=np.mean(lipstdmicelleradius)
meanpepcm=np.mean(meanpepcm,axis=0)
meanrescm=np.mean(meanrescm,axis=0)
stdrescm=np.mean(stdrescm,axis=0)
print(meanrescm)
xres=np.arange(1,len(meanrescm)+1)
pepcolor='#06d6a0'
miccolor='#ffd166'
plt.plot(xres,meanrescm,color=pepcolor,label='peptide')
plt.fill_between(xres, meanrescm-stdrescm, meanrescm+stdrescm,color=pepcolor, alpha=0.5)
plt.axhline(y=meanpepcm, linestyle='--',color=pepcolor)
plt.axhline(y=meanmicelleradius, color=miccolor, linestyle='--',label='micelle')
plt.fill_between(xres, meanmicelleradius-lipstdmicelleradius, meanmicelleradius+lipstdmicelleradius,color=miccolor, alpha=0.5)
plt.legend()

plt.title('ifp-lmpg')
plt.xticks(xres,resnames,fontsize=8, rotation=60)
plt.ylabel("Distance (Å)")
plt.ylim(0,35)
plt.xlim(1,len(xres))
plt.show()
