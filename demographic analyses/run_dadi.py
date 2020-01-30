import numpy as np
from numpy import array
import pylab
import dadi
import matplotlib.pyplot as plt



#function to convert generic sfs to dadi formatted text file
def fsToDadiFormat(arr):
    return str(len(arr)) + "\n" + " ".join([str(i) for i in arr])



file_BFAL = '/group/sbs_ssin/stella/albatross/05_reseq/initial_pipeline/06f_ANGSD/02_safQ/BFAL_autosome_1DSFS.sfs'
file_LAAL = '/group/sbs_ssin/stella/albatross/05_reseq/initial_pipeline/06f_ANGSD/02_safQ/LAAL_autosome_1DSFS.sfs'

fsArray_BFAL = np.loadtxt(file_BFAL)
fsArray_LAAL = np.loadtxt(file_LAAL)


#total number of sites is necessary for calculating scaling parameters
L_BFAL = int(round(np.sum(fsArray_BFAL),0))
L_LAAL = int(round(np.sum(fsArray_LAAL),0))

with open(file_BFAL.rstrip(".sfs") + ".dadi.sfs", "w") as ff: ff.write(fsToDadiFormat(fsArray_BFAL))
with open(file_LAAL.rstrip(".sfs") + ".dadi.sfs", "w") as ff: ff.write(fsToDadiFormat(fsArray_LAAL))


#load sfs as dadi object
data_BFAL = dadi.Spectrum.from_file(file_BFAL.rstrip(".sfs") + ".dadi.sfs")
data_LAAL = dadi.Spectrum.from_file(file_LAAL.rstrip(".sfs") + ".dadi.sfs")

data_BFAL.sample_sizes
data_LAAL.sample_sizes


#plot SFS for each population
plt.rcParams['figure.figsize'] = [10,12]
plt.subplot(2,1,1)
plotSpectrumBar(data_BFAL[~data_BFAL.mask], width=0.8, col = ["black"])
plt.subplot(2,1,2)
plotSpectrumBar(data_LAAL[~data_LAAL.mask], width=0.8, col = ["black"])


#estimate pi and Tajima's D from the SFS

print ("pi BFAL:", round(data_BFAL.pi()/L_BFAL,4))
print ("pi LAAL:", round(data_LAAL.pi()/L_LAAL,4))

print ("Tajima's D BFAL:", data_BFAL.Tajima_D())
print ("Tajima's D LAAL:", data_LAAL.Tajima_D())

