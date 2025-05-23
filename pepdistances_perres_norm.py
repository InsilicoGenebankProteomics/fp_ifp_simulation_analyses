import matplotlib.pyplot as plt
import numpy as np
import sys

nfiles=len(sys.argv)-1
meannormpepcm=[]
meannormrescm=[]
stdnormrescm=[]
lipstdmicelleradius=[]
for j in range(1,int(nfiles)+1):
	cmdistancesfile=str(sys.argv[j])
	print(cmdistancesfile)
	file = open(cmdistancesfile)
	lines = file.readlines()
	file.close()

	nlines = len(lines)
	nvars = len(lines[2].strip().split())
	normrescm = np.zeros((nlines-2,nvars-4))
	rescm = np.zeros(nvars-4)
	print(nlines,nvars)
	for i in range(2,nlines):
		l=lines[i].strip().split()
		rescm[:] = l[4:]
		normrescm[i-2,:]=rescm/float(l[1])
		lipstdmicelleradius.append(float(l[2])/float(l[1]))
	normpepcm=np.mean(normrescm,axis=1)
	meannormpepcm.append(np.mean(normpepcm))
	meannormrescm.append(np.mean(normrescm,axis=0))
	stdnormrescm.append(np.std(normrescm,axis=0))
meannormpepcm=np.mean(meannormpepcm,axis=0)
meannormrescm=np.mean(meannormrescm,axis=0)
stdnormrescm=np.mean(stdnormrescm,axis=0)
lipstdmicelleradius=np.mean(lipstdmicelleradius)
print(meannormrescm)
xres=np.arange(1,len(meannormrescm)+1)
pepcolor='#06d6a0'
miccolor='#ffd166'

plt.plot(xres,meannormrescm,color=pepcolor,label='peptide')
plt.fill_between(xres, meannormrescm-stdnormrescm, meannormrescm+stdnormrescm,color=pepcolor, alpha=0.5)
plt.axhline(y=meannormpepcm, linestyle='--',color=pepcolor)
plt.axhline(y=1, color=miccolor, linestyle='--',label='micelle')
plt.fill_between(xres, 1-lipstdmicelleradius, 1,color=miccolor, alpha=0.5)
plt.legend()

plt.title('fp')
#resnames=['1G','2A','3A','4L','5Q','6I','7P','8F','9A','10M','11Q','12M','13A','14Y','15R','16F']
resnames=['1M','2W','3K','4T','5P','6T','7L','8K','9Y','10F','11G','12G','13F','14N','15F','16S','17Q','18I','19L']
plt.xticks(xres,resnames,fontsize=8, rotation=60)
plt.ylabel("Normalized distance")
plt.ylim(0,1.2)
plt.xlim(1,len(xres))
plt.show()