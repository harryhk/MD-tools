# this code calculates the < \delta x(t) ^ 2 > between molecule i and molecule j 
# always use diffusion_thread.py without thread flag for single process. I keep this file for historical reason.

import xdrfile as xdrfile 
import numpy as np 
import sys , copy 
from  common.lnx_util import parseInput, print_help

_debug = 0 

class Diffusion(object):
    
    def __init__(self, inputP ):
        # class variable 
        # self.traj[i, j, k ]: i , ith frame 
        #                      j , jth atom 
        #                      k , 0 1 2 -- x, y, z 
        # self.molPair   : molPair to calcualte relative diffusion
        # self.timeStep  : timeStep between frames 
        # self.skip      : skip # of frames 
        # self.fout      : file handler of output file 
        # self.nFrames   : total # of frames in traj 
        # self.deltaF    : [] of delta frames for diffusion 
       
        trajFn = inputP['-f']
        molFn  = inputP['-n']


        # read traj 
        tmpx = []
        # tmpbox = []
        tmpT1 = -1 
        tmpT2 = -1 
        
        for f in xdrfile.xdrfile(trajFn):
            tmpx.append( copy.deepcopy(f.x)  )
            #tmpbox.append( copy.deepcopy(f.box) )
            
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
                if  f.time > 80.0 :
                    break 


        self.traj = np.array( tmpx )
        self.nFrames = np.size(self.traj , 0 ) 

        #self.box  = np.array( tmpbox )

        self.fout = open( inputP['-o'], 'w' ) 
        if inputP.has_key('-skip'):
            self.skip =  int( inputP['-skip'] ) 
        else:
            self.skip = 1 
        
        # nFrames / 10 will only provide 10 data points. 
        self.deltaF = range(1, self.nFrames/ 10, self.skip ) 
        
        if _debug:
            #print self.box
            print self.traj[:, 0, :]

        
        # read mol index 
        # format i ,j start from index 0 
        self.molPair = []
        for i in open(molFn):
            tmp1 , tmp2 = map(int ,  i.split() ) 
            self.molPair.append( (tmp1 -1, tmp2 -1) )

        if _debug:
            print self.molPair
            print self.timeStep
    
        

    def msd(self, moli, molj ):
        # r( t  ) should be centered by molecule i 
        # d( \delta t ) =  mean( (   r(t + \delta t  ) - r(t)    ) ^2 )

        # only x , y dimension is used 

        trajMoli =  self.traj[:, moli, 0:2 ] 
        trajMolj =  self.traj[:, molj, 0:2 ] 
        
        if _debug:
            print trajMoli
            print trajMolj

        deltaR  = trajMolj - trajMolj[0, :] - ( trajMoli  - trajMoli[0, : ] )  

        if _debug:
            print deltaR 
        
        msdList = [] 

        
        for deltaf in self.deltaF:
            tmp = 0  
            for i in range( self.nFrames - deltaf   ):
                nexti = i + deltaf 
                dr = deltaR[nexti, :] - deltaR[i, :] 
                tmp += np.dot( dr , dr )  
            

            msdList.append( tmp / (self.nFrames - deltaf )    ) 

        
        return msdList

    def allPair_msd(self):
        # get all pair msd from self.molPair
        tmp = []
        tmp.append( [ i * self.timeStep  for i in self.deltaF ] )

        for i, j in self.molPair:
            tmp.append( self.msd(i, j  ) )

        # print header 
        headerstr = '%15s' %('# time') + ''.join( [ '%15s' %   str(i)+'-'+str(j)   for i, j in self.molPair  ] ) 
        print >> self.fout , headerstr

        # transpose list 

        tmp = map( list, zip(*tmp) ) 
        for i in tmp:
            for j in i : 
                self.fout.write( '%15.5f' % j ) 

            self.fout.write('\n')


            
        

if __name__ == '__main__':
    
    inputP = parseInput(sys.argv[1:])

    paraOpt = '-f -n -o -h -skip'.split(' ')

    helpdoc = ''\
              'Usage prog.py  -f trajNojump.xtc ; traj of nojump from g_traj  \n'\
              '               -n mol.ndx  ; relative diffusion between mol i and mol j , index start from 1 \n'\
              '               -skip  1[defualt]    ; skip frames in delta T \n'\
              '               -o out.xvg  ; \n'
    
    print_help(inputP, paraOpt , helpdoc )

    prob = Diffusion( inputP )

    prob.allPair_msd()
