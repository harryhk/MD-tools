# this code calculate the all pair wise distance. 

import xdrfile as xdrfile 
import numpy as np 
import sys , copy 
from  common.lnx_util import parseInput, print_help

_debug = 0 

class Survival(object):
    
    def __init__(self, inputP ):
        # class variable 
        # self.traj[i, j, k ]: i , ith frame 
        #                      j , jth atom 
        #                      k , 0 1 2 -- x, y, z 
        # self.box[i, j] : i, ith frame 
        #                : j, box x, y ,z dimention 
        # self.timeStep  : timeStep between frames
        # self.molPair   : initial pairs  
        # self.nFrames   : total # of frames in traj 
        # self.inputP    : inputP ; keep it will be convinient 

        trajFn = inputP['-f']
        molFn  = inputP['-n']

        self.inputP  = inputP 


        # read traj 
        tmpx = []
        tmpbox = []
        tmpT1 = -1 
        tmpT2 = -1 
        

        for f in xdrfile.xdrfile(trajFn):
            tmpx.append( copy.deepcopy(f.x)  )
            tmpbox.append( np.diag( copy.deepcopy(f.box) ) )
            
            if tmpT1 == -1 and tmpT2 == -1:
                tmpT1 = f.time 
            elif tmpT1 != -1 and tmpT2 == -1:
                tmpT2 = f.time 
                self.timeStep = tmpT2 - tmpT1
            else:
                tmpT1 = tmpT2 
                tmpT2 = f.time 

                if self.timeStep != tmpT2 - tmpT1 :
                    print >> sys.stderr, 'Time step changed ! '
                    sys.exit(1)

            if _debug: 
                if  f.time > 20.0 :
                    break 


        self.traj = np.array( tmpx )
        self.nFrames = np.size(self.traj , 0 ) 

        self.box  = np.array( tmpbox )

        if _debug:
            pass 
            #print self.box
            #print self.traj

        
        # read mol index 
        # format i ,j start from index 0 
        self.molPair = []
        for i in open(molFn):
            tmp1 , tmp2 = map(int ,  i.split() ) 
            self.molPair.append( (tmp1 -1, tmp2 -1) )

        if _debug:
            #print self.molPair
            print self.timeStep
        
    
    def dist( self, moli, molj, startIdx ):
        
        coorMoli = self.traj[ startIdx , moli , 0:2 ]
        coorMolj = self.traj[ startIdx , molj , 0:2 ] 

        box = self.box[startIdx][0:2]
        
        tmp = coorMoli - coorMolj 

        tmp = tmp - np.round(  tmp / box )  * box 

        dist = np.sqrt( np.dot( tmp, tmp ) ) 

        return dist


    def singleRate(self ):

        # determine initially mol-mol pairs within cutoff
        # startIdx * self.timeStep is start time 

        data = np.zeros( ( self.nFrames, len( self.molPair )    ) )
        
        for timeIdx in range(self.nFrames):
            
            for idx, (i, j) in enumerate( self.molPair ) :
                 
                data[timeIdx , idx ] = self.dist( i, j, timeIdx ) 
                

            if self.inputP.has_key('-v'):
                sys.stderr.write('Progress: %3.2f\r' % (timeIdx / float( self.nFrames )  ) )
                sys.stderr.flush() 




        # save data 
        fout = open( self.inputP['-o'] , 'w')
        np.save( fout , data ) 
        fout.close() 


    def solv(self):
        self.singleRate() 

if __name__ == '__main__':
    
    inputP = parseInput(sys.argv[1:])

    paraOpt = '-f -n -h  -o -v '.split(' ')

    helpdoc = ''\
              'Usage prog.py  -f trajNojump.xtc ; traj of nojump from g_traj  \n'\
              '               -n mol.ndx  ; pair wise distance between mol i and mol j , index start from 1 \n'\
              '               -v               ; output progress \n'\
              '               -o   bit.npy     ; numpy bit array to all pair cutoff\n'
    
    print_help(inputP, paraOpt , helpdoc )

    prob = Survival( inputP )

    prob.solv()
