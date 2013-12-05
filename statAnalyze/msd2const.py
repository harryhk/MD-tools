# this code take msd from diffusion_thread.py and extract diffusion constant 

import numpy as np 
import sys
from common.lnx_util import parseInput, print_help 


if __name__ == '__main__':
    

    inputP = parseInput(sys.argv[1:])

    paraOpt = '-f -o -h'.split()

    helpdoc = ''\
              'Usage prog.py  -f *.xvg  ; all xvg files that contains msd\n'\
              '               -o  out.xvg  ; output file \n'

    print_help(inputP , paraOpt, helpdoc )

    result= {}  # {'0-1' : diffusion constant }
    
    msdfns = inputP['-f']
    if type(msdfns) != list :
        msdfns = [ msdfns ]
        

    for msdfn in msdfns: 
        
        tmpdata = np.loadtxt(msdfn)

        header = open(msdfn).readline().strip().split()[2:]

        if len(header) != np.size(tmpdata, 1) -1  :
            print header 
            print len(header) , np.size( tmpdata, 1 ) 
            print >> sys.stderr, "data not match "
            sys.exit(1)
        
        for i in range( len( header) ) :     
            p = np.polyfit(tmpdata[:,0], tmpdata[:,i+1], 1)

            result[ header[i] ] = p 


    fout = open(inputP['-o'], 'w')
    for i in sorted( result.keys() , key=lambda x : map( int, x.split('-') ) ) :
        print >> fout , '%15s%15.5g%15.5g' % (i, result[i][0], result[i][1])

     
        
