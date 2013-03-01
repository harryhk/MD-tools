#!/usr/bin/env python 

import sys, pdb
import numpy as np

def xvg( fn  ):
	f=open(fn, 'r')
	return np.array( [  i.strip().split()  for i in f if i[0] !='@' and i[0] != '#'] , dtype=float)


if len(sys.argv) != 4:
	print("Usage! prog.py dens.xvg  water_dens(am) [am/nm]")
	exit(1)

#water_s=float(sys.argv[3])
#water_e=float(sys.argv[4])

data=xvg(sys.argv[1])

#data_w = np.array(  [ i[2] for i in data if i[0] > water_s and i[0] < water_e ]     )
water_dens = float( sys.argv[2] )

# convert nm to am 
if sys.argv[3] =='nm':
	data = np.array( [ [ i[0]*10, i[1] /1000 ] for i in data   ] )

# substract bulk water density 
data[:,1] = data[:,1] - water_dens

#pdb.set_trace()

# calculate dz , since dz is even distributed 
dz =  data[1:, 0] - data[0:-1, 0]
dz = dz.mean() 
z  = ( data[1:,0] + data[0:-1, 0] ) / 2.0

# calcuate \rho(z)
av = ( data[1:,1] + data[0:-1, 1] ) / 2.0 

q = np.linspace(0, 1 , num = 1000)

cos_q =   np.outer( z, q  )
cos_q = np.cos(cos_q)

temp_sum = av * cos_q.transpose()

temp_sum = temp_sum.transpose()

sum = np.sum( temp_sum, 0) * dz 

for i , j in zip( q, sum  ):
	print "%10.3f%10.3f" % ( i, np.abs(j))

#for q in np.linspace(0,1, num=1000):
#	sum=0;
#
#	for i in range(np.size(data,0) -1 ) :
#		av= ( data[i, 1] + data[i+1, 1] ) / 2
#		dz= data[i+1,0] - data[i,0]
#		z= (data[i+1,0] + data[i,0] ) /2
#		sum+= av * np.cos( q * z ) * dz
#	
#
#
#	print "%10.3f%10.3f" % ( q, 2*np.abs(sum)  )  


		


