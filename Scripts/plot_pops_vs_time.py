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
recs = np.fromfile(f, dtype=rectype, offset=8)

timestep = []
for i in range(1, len(recs['time'])):
	timestep.append( recs['time'][i] - recs['time'][i-1])
timestep = np.array(timestep)

plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,0]), label="level id = 1")
plt.plot(np.log10(recs['time'][1:]), np.log10(recs['pops'][1:,2]), label="level id = 3")

plt.legend()
plt.show()
