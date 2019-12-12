#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

f = open('pops_vs_time_6.bin')
num_of_levs = np.fromfile(f, dtype=np.dtype('u8'), count=1)
num_of_levs = int(num_of_levs[0])
rectype = np.dtype( [ ('time', 'd'), ('pops', 'd', num_of_levs) ] )
recs0 = np.fromfile(f, dtype=rectype)

timestep0 = []
for i in range(1, len(recs0['time'])):
	timestep0.append( recs0['time'][i] - recs0['time'][i-1])
timestep0 = np.array(timestep0)

#plt.plot(recs['time'], recs['pops'][:,766])
plt.plot(np.log10(recs0['time'][1:]), np.log10(recs0['pops'][1:,1827]), label="old")
plt.plot(np.log10(recs0['time'][1:]), np.log10(recs0['pops'][1:,624]), label="old")
#plt.plot(np.log10(recs['time'][1:]), np.log10(timestep), 'o')
#for i in range(1, num_of_levs-1):
#for i in range(num_of_levs-1, 1, -1):
	#plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,i]/recs['pops'][1:,0]) + (num_of_levs-i)*0.01)
	#plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,i]/recs['pops'][1:,0]) + np.log10(recs['pops'][1:,i]/recs['pops'][1:,0])[0])
	#plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,i]), 'o')

f = open('pops_vs_time_NR.bin')
num_of_levs = np.fromfile(f, dtype=np.dtype('u8'), count=1)
num_of_levs = int(num_of_levs[0])
rectype = np.dtype( [ ('time', 'd'), ('pops', 'd', num_of_levs) ] )
recs = np.fromfile(f, dtype=rectype)

timestep = []
for i in range(1, len(recs['time'])):
	timestep.append( recs['time'][i] - recs['time'][i-1])
timestep = np.array(timestep)

#plt.plot(recs['time'], recs['pops'][:,766])
plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,1827]))
plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,624]))
#plt.plot(np.log10(recs['time'][1:]), np.log10(timestep), 'o')
#for i in range(1, num_of_levs-1):
#for i in range(num_of_levs-1, 1, -1):
	#plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,i]/recs['pops'][1:,0]) + (num_of_levs-i)*0.01)
	#plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,i]/recs['pops'][1:,0]) + np.log10(recs['pops'][1:,i]/recs['pops'][1:,0])[0])
	#plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,i]), 'o')

plt.legend()
plt.show()
