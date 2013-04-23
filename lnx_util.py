# this file provides basic utilities for other analyse tools. 
import re
import sys


def parseInput(l):
    ''' l = sys.argv[1:]. this will parse the input and return a dictionary.
    '''
    
    inputP = {}
    
    flag = re.compile('-\D+')
    indx = [ i for (i, n ) in enumerate(l) if flag.match(n)   ]
    indx.append(len(l))
    
    for i in range( len(indx) -1 ):
        inputP[l[indx[i]]] = l[ indx[i]+1: indx[i+1] ]
        if len(inputP[l[indx[i]] ] ) ==1 :
            inputP[ l[indx[i]]] = inputP[l[indx[i]]][0]

    return inputP


def print_help(inputP, paraOpt, helpdoc):
    ''' inputP: dictionary parsed from parseInput
        paraOpt: legit keys in inputP 
        helpdoc: help string, trigerred by -h  
    '''
    
    for i in inputP:
        if i not in paraOpt or i == '-h':
            sys.stderr.write(helpdoc)
            sys.exit(1)

