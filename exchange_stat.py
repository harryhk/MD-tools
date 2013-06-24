#!/bin/env python 

# calculate the replica exchange statistics 

import sys, re
from lnx_util import parseInput , print_help

inputP= parseInput(sys.argv[1:])
paraOpt = ' -fin  -n  -h '.split()

helpdoc=\
'''
Usage! ./prog.py 
       -fin exchange.dat
'''

print_help(inputP, paraOpt, helpdoc)

fin = [ i.strip() for i in open(inputP['-fin']) ]
fin = [ re.split('\d', i)[1:-1] for i in fin   ]

nrepl = len(fin[0] )
n_attempt = len(fin)

stat = [ 0 for i in range(nrepl) ]
for i in fin:
    for jc, j in enumerate(i):
        if 'x' in j :
            stat[jc] +=1

# print statistics
print '%5d     '* (nrepl +1 ) % tuple(range(1, nrepl +2))
print '     %5.2f' * nrepl  % tuple( [ float(i) / n_attempt * 2.0 for i in stat  ]   )



