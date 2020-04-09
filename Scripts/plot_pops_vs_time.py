#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

#f = open(sys.argv[1])
f = open('pops_vs_time.bin')
num_of_levs = np.fromfile(f, dtype=np.dtype('u8'), count=1)
num_of_levs = int(num_of_levs[0])
rectype = np.dtype( [ ('time', 'd'), ('pops', 'd', num_of_levs) ] )
recs = np.fromfile(f, dtype=rectype)

timestep = []
for i in range(1, len(recs['time'])):
	timestep.append( recs['time'][i] - recs['time'][i-1])
timestep = np.array(timestep)

plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,1827]), label="level id = 1828")
plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,624]), label="level id = 625")

plt.legend()
plt.show()
