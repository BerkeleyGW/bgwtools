#!/usr/bin/env python

from __future__ import division
from numpy import *

def get_data_from_eqp(fname, column=1):
	X  = []
	Y  = []
	Z  = []
	Vl_LDA = []
	Vl_GW  = []
	#read from eqp file
	f=open(fname)
	item_num=-1
	for line in f.readlines():
		items=line.split()	
		if (item_num<0):
			kpts=map(float,items[0:3])
			#Kpts+=kpts
			X.append(kpts[0])
			Y.append(kpts[1])
			Z.append(kpts[2])
			items_cnt=int(items[3])
			Vl_LDA.append({})
			Vl_GW.append({})
			item_num=0
		else:
			band_idx=int(items[1])
			energ_LDA=float(items[2])
			energ_GW =float(items[3])
			Vl_LDA[-1][band_idx]=energ_LDA
			Vl_GW[-1][band_idx] =energ_GW
			item_num+=1
			if (item_num==items_cnt):
				item_num=-1

	bands=[v.keys() for v in Vl_LDA]
	bands=array(sorted(unique(array(bands))))
	#print bands
	V_LDA = ma.masked_all( (len(X),len(bands)) )
	V_GW  = ma.masked_all( (len(X),len(bands)) )
	i=0
	for v_ in Vl_LDA:
		indices = [where(bands==idx)[0][0] for idx in v_.keys()]
		V_LDA[i,:][indices] = v_.values()
		V_LDA.mask[i,:][indices] = False
		i+=1
	i=0
	for v_ in Vl_GW:
		indices = [where(bands==idx)[0][0] for idx in v_.keys()]
		V_GW[i,:][indices] = v_.values()
		V_GW.mask[i,:][indices] = False
		i+=1

	return bands,array(X),array(Y),array(Z), V_LDA, V_GW

