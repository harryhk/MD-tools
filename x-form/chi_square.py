#!/user/bin/env python 
# this code calculate the chi square value between experimental form factor and simulation one
# make sure that your experimental file format :    
#    q , f(q) , delta f(q)  
# simulation form factor format 
#    q , f(q) 
# the code linear fit the simulation f(q) to experimental data points 

import sys
import numpy as np 
from common.lnx_util import parseInput , print_help, xvg 

inputP = parseInput(sys.argv[1:])
paraOpt = ' -fexp  -fsim   '.split()

helpdoc=\
'''
Usage! ./prog.py
       -fexp  file1 ; input file for experimental form factor
       -fsim  file2 ; input file for simulation form factor 
'''

print_help(inputP, paraOpt, helpdoc)

fexp = xvg(inputP['-fexp'])
fsim = xvg(inputP['-fsim'])

chi = 0 
i=0
j=1
for q, fq, fqdelta in fexp:
    # use linear interpolation to find simulated value 
    while not ( fsim[i,0] < q and fsim[j,0] >= q ):
        i+=1
        j+=1

    #print q, fsim[i,0], fsim[j,0]
    fqsim = ((fsim[j,0]-q) * fsim[i,1]  + ( q-fsim[i,0]) * fsim[j,1] )/(fsim[j,0] - fsim[i,0])
    
    tmp = ( (np.abs(fq) - np.abs(fqsim) ) / fqdelta )** 2
    #print  q, fq, fqsim, fqdelta , tmp
    chi += tmp

chi /= len(fexp) -1
chi = np.sqrt( chi) 
print chi

        


