#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

class data_class:
	def _init(self):
		self._num_of_species = 0
		self._num_of_trans = []
		self._recs = []

	def read_data(self, filename):
		with open(filename) as f:
			self._num_of_species = np.fromfile(f, dtype=np.int64, count=1)[0]
			self._num_of_trans = np.fromfile(f, dtype=np.int64, count=self._num_of_species)
			tot_num_of_trans = np.sum(self._num_of_trans, dtype=np.int64)
			rectype = np.dtype( [ ('time', np.double), ('pops', np.double, tot_num_of_trans) ] )
			#self._recs = np.memmap(f, dtype=rectype, mode='r')
			self._recs = np.fromfile(f, dtype=rectype)
	
	def __init__(self, filename=None):
		self._init()
		if filename is not None:
			self.read_data(filename)

	def get_time(self):
		return self._recs['time'][1:]

	def get_data(self, ispec, itrans):
		trans_to_begin = np.sum(self._num_of_trans[0:ispec], dtype=np.int64)
		return self._recs['pops'][1:, trans_to_begin + itrans]

def find_gradient(t, f):
	res = []
	time = []
	for i in range(1, len(f)):
		df_dt = (f[i] - f[i-1]) / (t[i] - t[i-1])
		res.append(df_dt)
		res.append(df_dt)
		time.append(t[i-1])
		time.append(t[i]*0.999)

	return np.array(time), np.array(res)


#recs = data_class(sys.argv[1])
recs = data_class('pops_vs_time.bin')

plt.plot(recs.get_time(), np.log10(recs.get_data(0, 0)), label="p-H2CO level id = 1")
plt.plot(recs.get_time(), np.log10(recs.get_data(1, 2)), label="o-H2CO level id = 3")

plt.legend()
plt.show()
