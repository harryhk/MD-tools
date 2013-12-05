#!/usr/bin/env python 
import sys
import numpy as np

# read from stdin and calcualte statistics 

d = np.array( [ i.strip() for i in sys.stdin.readlines()  ]  , dtype=float)

print d.mean() , d.std()/ np.sqrt( np.size(d, 0 )   ) 



