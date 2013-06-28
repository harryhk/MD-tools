#!/usr/bin/env python 

import sys, pdb
import numpy as np 

#data= np.array( open('t').readlines() , dtype=float) 

data1 = np.random.normal(size=100000)

data= np.zeros_like( data1 )
data[0] = data1[0]
rho = 1/1.01

for i in range(1, len(data1)):
	data[i] = rho * data[i-1] + np.sqrt( 1- rho*rho ) * data1[i]

n = len(data) 
n_each_block=1

while n_each_block <= n/100:
	
	n_blocks = n / n_each_block; 
	
	data_new = data[0:n_blocks * n_each_block].reshape( n_blocks, n_each_block  )
	data_new = np.mean(data_new, 1)

	
	c0 = np.var(data_new)
	print "%10.3f%10.3f" % ( np.sqrt( c0 / (n_blocks-1) )  , np.sqrt(c0/2.0) /(n_blocks-1) ) 
	n_each_block = n_each_block + 10


