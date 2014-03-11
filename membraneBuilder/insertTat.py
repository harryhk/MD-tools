# this script used to insert Tat into desired depth 

import numpy as np 
import sys, copy
import common.lnx_util as lu 

class Index(object):
    
    def __init__(self, lines)
        



class Gro(object):
    
    def __init__(self, lines):
        
        self.numAtom = len(lines)
        self.resid = [  int(i[:5]) for i in lines ]
        self.resname =  [ i[5:10].strip() for i in lines  ]
        self.atomName = [ i[10:15].strip() for i in lines ]
        self.index = np.array( [ i[15:20] for i in lines] ,dtype=int )
        self.coor = np.array([ i[20:44].split() for i in lines ], dtype= float)

    def selectByIndex(self, index):
        
        # select coor by index, 
        # notice here we are not doing copy 
        return self.coor[index-1]



