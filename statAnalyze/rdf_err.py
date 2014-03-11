#!/usr/bin/env python 

import sys
from common.lnx_util import parseInput , print_help
import numpy as np 

# rdf range [0, 3.78] ; please run rdf_range.sh to check 

_debug = 0

def readfile(filen, rmin, rmax):
    
    tmp =  [ map( float , i.strip().split() ) for i in open(filen)  if i[0] != '@' and i[0] != '#' ] 
    tmp = np.array( [  i for i in tmp if i[0]<=rmax and i[0] >= rmin   ] ) 
    return tmp

def absErr(l1, l2):
    
    if _debug:
        print 'l1\n', l1[-10:-1]
        print 'l2\n', l2[-10:-1]
        print abs(l1-l2)[-10:-1]
    return sum(abs( l1-l2) ) / len(l1)

def weightErr(l1,l2):
    
    tmp = ( l1 + l2 ) /2.0 
    c= 0 
    for i in range( len(tmp) ):
        if tmp[i] == 0 :
            continue
        else:
            c += abs(l1[i] - l2[i])/ tmp[i]

    return c / len(l1) 
    


inputP = parseInput( sys.argv[1:] ) 
paraOpt = '-f  -range  -abs  -h -r '.split(' ')

helpdoc = 'Usage ./prog.py  -f rdfu.xvg rdfl.xvg  ; rdf upper and lower files \n'\
          '                 -range 0 3.78         ; only calculate rdf among range \n'\
          '                 -abs                  ; if present , use absolute error of each point, if not, use weighted\n'\
          '                 -r                    ; if present, g(r) * r '

print_help(inputP, paraOpt, helpdoc)

rmin = float(inputP['-range'][0])
rmax = float(inputP['-range'][1])

finu = readfile(inputP['-f'][0], rmin, rmax) 
finl = readfile(inputP['-f'][1], rmin, rmax)

if _debug:
    print "finu\n", finu[-10:-1]
    print "finl\n", finl[-10:-1]

# check dimension of finu and finl


if np.size(finu, 0) != np.size( finl , 0 ) :
    print "# of rows not match "
    sys.exit(1)

if inputP.has_key( '-r' ) :
    finu[:,1] *= finu[:,0] * finu[:,1]
    finl[:,1] *= finl[:,0] * finl[:,1]

if inputP.has_key('-abs'):
    print absErr( finu[:,1], finl[:,1]  )
else:
    print weightErr( finu[:,1], finl[:,1]  )
