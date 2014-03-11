#!/usr/bin/env python 
# this script will choose a frame that has closest area per lipid as we specified. 
# It also considers the frequency that we save xtc files 


import sys
from common.lnx_util import parseInput, print_help, xvg

if __name__ == '__main__':
    
    inputP = parseInput( sys.argv[1:]  )  

    paraOpt = '-f  -area  -n  -mod -h'.split(' ')
    helpdoc =  'Usage ./prog.py -f box.xvg  ; box information from g_energy \n'\
               '-area 68.0    ;  find most close area per lipid to 68 a^2 \n'\
               '-n  64        ;  number of lipids per leaflet \n'\
               '-mod  10      ;  the frame should also be divided by 10 ps, since we are saving by every 10 ps\n'


    print_help(inputP, paraOpt, helpdoc)


    fin = xvg(inputP['-f']) 
    target = float(inputP['-area'])
    mod = int(inputP['-mod'])
    nlipids = int(inputP['-n'])

    apl = []
    for i in fin:
        apl.append( (i[0], abs(i[1]*i[2]/nlipids*100 - target  )  ) )

    
    apl.sort(key=lambda x:x[1]  ) 
    
    for i , data in enumerate(apl):
        if data[0] % mod == 0:
            print "%d ps, %f" %  (data[0] , data[1])
            break


