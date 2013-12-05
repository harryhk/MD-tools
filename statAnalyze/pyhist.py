#!/usr/bin/env

# simple histogram code ; 

from common.lnx_util import parseInput, print_help
import sys 
import numpy as np 

if __name__ == '__main__':
    
    inputP = parseInput(sys.argv[1:])

    paraOpt = '-xrange -bins -density'.split() 
    
    helpdoc =\
'''
Usage prog.py : programe read from stdin and write to stdout 
    -xrange  xmin, xmax ; if not present, use numpy default 
    -bins  ; default is 10 
    -density ; if set, output probability density. 
'''
    
    print_help(inputP, paraOpt , helpdoc) 

    r = None
    if inputP.has_key("-xrange"):
        r = tuple( inputP['-xrange']  )

    
    bins= 10 
    if inputP.has_key("-bins"):
        bins = int( inputP['-bins']) 

    dens = False
    if inputP.has_key("-density"):
        dens = True 

    data = [ i.strip() for i in sys.stdin.readlines() ] 
    data = np.array( data, dtype = float) 

    hist, bin_edges = np.histogram( data, bins, range= r, density = dens ) 

    
    for i in range(len( hist ) ):
        print "%15.5g%15.5g" % (( bin_edges[i] + bin_edges[i+1] ) / 2.0,  hist[i]  )
        
