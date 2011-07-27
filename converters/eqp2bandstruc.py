#!/usr/bin/env python

import sys
from numpy import *

fname = sys.argv[1]
f=open(fname)

item_num=-1
for line in f.readlines():
	items=line.split()	
	if (item_num<0):
		kpts=tuple(map(float,items[0:3]))
		kpts_str=('  %10.5f'*3)%kpts
		items_cnt=int(items[3])
		kpt_data=[]
		item_num=0
	else:
		p_str='%8d'%(int(items[1]))
		p_str+=kpts_str
		#kpt_data += [int(items[1])]+map(float,items[2:4])
		energies=map(float,items[2:4])
		energies=tuple(energies+[energies[1]-energies[0]])
		p_str+=('  %10.5f'*3)%energies
		print p_str
		item_num+=1
		if (item_num==items_cnt):
			item_num=-1

