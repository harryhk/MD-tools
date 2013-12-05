#!/usr/bin/env python 

# calcualte surface tension 

import sys
from common.lnx_util import parseInput , print_help , xvg 

import numpy as np 


inputP = parseInput(sys.argv[1:])
paraOpt = "-f -col -n  -h ".split(' ')

helpdoc = 'Usage ./prog.py  -f file.xvg  ; input \n'\
          '                 -col px py pz lz ; col indx start from 0 \n'\
          '                 -n  2  ; number of surfaces \n'\

print_help(inputP , paraOpt , helpdoc )

data = xvg(inputP['-f'])

colpx, colpy, colpz , collz  = map( int, inputP['-col'] ) 

nsurfs = int( inputP['-n'] )

surfT =  ( data[:, colpz] - ( data[:, colpx] + data[:, colpy] ) / 2.0 ) * data[:, collz]

surfT = surfT / nsurfs ; 

for i in surfT : 
    print "%15.5f" % i 



