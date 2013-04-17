#!/usr/bin/env python 

import sys, pdb
import numpy as np
from  lnx_util import parseInput, print_help

def xvg( fn  ):
    f=open(fn, 'r')
    return np.array( [  i.strip().split()  for i in f if i[0] !='@' and i[0] != '#'] , dtype='float')

inputP = parseInput(sys.argv[1:])
paraOpt = ['-f', '-nm2am' , '-wdens', '-q' , '-of' , '-oq', '-h']

helpdoc="Usage!   ./prog.py\n"\
        "         -f input.edens  \n"\
        "         -nm2am   ; convert unit from nm to am \n"\
        "         -wdens  0.323  ; bulk water density. This has to be done manually.\n"\
        "         -q  0.0 0.25 0.5 ... ; select multiple q value to look at each z contriubtion \n"\
        "         -of  out.xff   ; filename for output xray form factor\n"\
        "         -oq  out.qz    ; filename, contribution from each z at different q\n"\
        "         -h             ; print help document\n"

print_help(inputP, paraOpt, helpdoc)

#if len(sys.argv) != 4:
#    print("Usage! prog.py dens.xvg  water_dens(am) [am/nm]")
#    exit(1)

#water_s=float(sys.argv[3])
#water_e=float(sys.argv[4])

data=xvg(inputP['-f'])

#data_w = np.array(  [ i[2] for i in data if i[0] > water_s and i[0] < water_e ]     )
water_dens = float(inputP['-wdens'])

# convert nm to am 
if inputP.has_key('-nm2am'):
    data = np.array( [ [ i[0]*10, i[1] /1000 ] for i in data   ] , dtype='float')
    water_dens /= 1000.0

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
dq = q[1]-q[0]

cos_q =   np.outer( z, q  )

cos_q = np.cos(cos_q)
# cos_q
# cos(z1q1), cos(z1q2), ...
# cos(z2q1), cos(z2q2), ...
# cos(z3q1), cos(z3q2), ...


temp_sum = av * cos_q.transpose()
# temp_sum
# av(z1)*cos(z1q1), av(z2)*cos(z2q1), av(z3)*cos(z3q1), ...
# av(z1)*cos(z1q2), av(z2)*cos(z2q2), av(z3)*cos(z3q2), ...
# av(z1)*cos(z1q3), av(z2)*cos(z2q3), av(z3)*cos(z3q3), ...


temp_sum = temp_sum.transpose()
# temp_sum
# av(z1)*cos(z1q1), av(z1)*cos(z1q2), av(z1)*cos(z1q3), ...
# av(z2)*cos(z2q1), av(z2)*cos(z2q2), av(z2)*cos(z2q3), ...
# av(z3)*cos(z3q1), av(z3)*cos(z3q2), av(z3)*cos(z3q3), ...

if inputP.has_key('-oq'):
    f1=open(inputP['-oq'], 'w')
    f1.write('%s\n' % ('#'+' '.join(sys.argv)) ) 
    
    # select q  to write 
    selectq = inputP['-q']
    if type(selectq) is not list:
        selectq = [ selectq]
    qindx = [  int(float(i) / dq)   for i in selectq ]

    for i, zi  in zip(temp_sum, z):
        qsForz= ''.join(  map( lambda x: '%15.6f' % x ,  [  i[ii] for ii in qindx  ] ) ) 
        f1.write('%15.6f%s\n' % (zi, qsForz) )
    f1.close()
        

if inputP.has_key('-of'):
    f2 = open(inputP['-of'],'w')
    f2.write('%s\n' % ( '#'+' '.join(sys.argv)))
    sum = np.sum( temp_sum, 0) * dz 
    for i , j in zip( q, sum  ):
        f2.write( "%10.3f%10.3f\n" % ( i, np.abs(j)) ) 

#for q in np.linspace(0,1, num=1000):
#    sum=0;
#
#    for i in range(np.size(data,0) -1 ) :
#        av= ( data[i, 1] + data[i+1, 1] ) / 2
#        dz= data[i+1,0] - data[i,0]
#        z= (data[i+1,0] + data[i,0] ) /2
#        sum+= av * np.cos( q * z ) * dz
#    
#
#
#    print "%10.3f%10.3f" % ( q, 2*np.abs(sum)  )  


        


