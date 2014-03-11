#!/user/bin/env python 

# this script caclualtes cdf from rdf 


import sys
from common.lnx_util import parseInput , print_help, xvg 
import numpy as np 
import scipy.integrate   


if __name__ == '__main__':
    
    inputP = parseInput(sys.argv[1:] ) 
    paraOpt = '-rdf -o -h   -d  '.split() 

    helpdoc = \
'''
Usage ./prog.py -rdf   ; rdf file 
                -o   ; output cdf ; if no present, sys.stdout 
                -d    ; dimension 
'''

    print_help( inputP, paraOpt, helpdoc ) 

    rdfdata = xvg( inputP['-rdf'] ) 

    r= rdfdata[:,0]
    g= rdfdata[:,1]


    if int( inputP['-d'] )  == 2 :
        cdf = scipy.integrate.cumtrapz( g * 2 * scipy.pi * r , r )
    elif int( inputP['-d'] == 3) :
        pass 
    else:
        print >> sys.stderr, "Dimension not supported"
        sys.exit(1) 
        
    
    if inputP.has_key('-o'):
        fout = open( inputP['-o'] , 'w' ) 
    else:
        fout = sys.stdout 

    for i , j in zip( r, cdf ) :
        print >> fout ,"%15.5f%15.5f" % (i, j) 
