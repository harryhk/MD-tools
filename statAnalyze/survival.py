# this code calculate the survival rate of mol-mol pairs. 
# it randomly choose nsample of starting point , and calcualte survial pairs vs time 

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
        # self.skip      : number of frames to skip in survial time  
        # self.cutoff    : cutoff 
        # self.nFrames   : total # of frames in traj 
        # self.inputP    : inputP ; keep it will be convinient 
        # self.nsample   : number of random start points draw from first half of the data   

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

        self.nsample = int( inputP['-nsample'] ) 
        
        self.skip = int( inputP['-skip'] )
        self.cutoff = float( inputP['-cutoff'] )
        
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
        
    
    def withinCutoff( self, moli, molj, startIdx ):
        
        coorMoli = self.traj[ startIdx , moli , 0:2 ]
        coorMolj = self.traj[ startIdx , molj , 0:2 ] 

        box = self.box[startIdx][0:2]
        
        tmp = coorMoli - coorMolj 

        tmp = tmp - np.round(  tmp / box )  * box 

        dist = np.sqrt( np.dot( tmp, tmp ) ) 

        if dist > self.cutoff:
            return False
        else:
            return True 



    def singleRate(self , startIdx ):

        # determine initially mol-mol pairs within cutoff
        # startIdx * self.timeStep is start time 

        pairs=[]
        sv = [] 
        
        for i, j in self.molPair:
            if self.withinCutoff(i, j , startIdx ):
                pairs.append( ( i ,j) ) 

        sv.append( len( pairs)  ) 

        if _debug:
            pass
            #print len(pairs)
            #print pairs[0]
            #print self.traj[startIdx][pairs[0][0]]
            #print self.traj[startIdx][pairs[0][1]]
            #print self.box[startIdx]
    
        
        #print "init \n" , pairs 
        # time from startIdx to self.nFrames with self.skip 

        # debug the survival 
        # debugPair = [(72, 93) , (83, 109), (124, 143) ]

        
        for timeIdx in range( startIdx + self.skip, self.nFrames, self.skip):
            
            c = 0 
            #tmppair = []
            for i, j in pairs:
                if self.withinCutoff( i, j, timeIdx ) :
                    c += 1 
                    #tmppair.append((i,j ))

            sv.append(c ) 

            #for i, j  in debugPair:
            #    sys.stdout.write( '%15.5f' % ( timeIdx * self.timeStep )   ) 
            #    if self.withinCutoff( i, j , timeIdx) :
            #        sys.stdout.write( '%15d' % 1  ) 
            #    else:
            #        sys.stdout.write( '%15d' % 0  ) 
            #sys.stdout.write( '\n' ) 
        
        # print output 
        #print 'end \n', tmppair 

        #sys.exit(1)
        #t = sv[0]
        #for i in sorted( sv.keys() ):
        #    print >> self.fout , '%15.5f%15d%15.5f' % (i*self.timeStep , sv[i] , sv[i] / float( t)  )


        return sv 
        
        

    def nSample(self):
        # get all pair msd from self.molPair
        
        #self.singleRate( 0) 

        nResults = [] 

        endIdx = self.nFrames / 2 
        
        # equally spaced start point from [0 , n/2 ] 
        c = 0 
        for startIdx  in map( int, np.linspace(0, endIdx , self.nsample ) ) : 
            
            nResults.append( self.singleRate(startIdx ) )
            c += 1
            if inputP.has_key('-v'):
                sys.stdout.write("Progress: %2.1f \r" %  ( float(c) / self.nsample )  )
                sys.stdout.flush() 


        return nResults 

    


    def solv(self):
        nResults = self.nSample() 
        
        # print nResults to file 
        if self.inputP.has_key('-o'):
            fout = open(inputP['-o'], 'w')
            for i in range(self.nsample):
                # print header 
                print >> fout , "# Sample %d" % i 
                t0 = nResults[i][0] 
                for j in range( len( nResults[i])  ):
                    fout.write('%15d%15d%15.5f\n' % (j * self.skip * self.timeStep , nResults[i][j] , nResults[i][j] / float( t0 ) ) )

                fout.write('\n\n')
            fout.close() 

        # get the length of each sample 
        nsampleLen = [ len(i)  for i in nResults ] 
        maxLen = max( nsampleLen ) 

        fout = open(inputP['-ostat'] , 'w')
        for i in range(maxLen ):
            tmp = []
            for j in range(self.nsample):
                if i < nsampleLen[j]:
                    tmp.append( nResults[j][i] / float( nResults[j][0] )  ) 

            tmp = np.array( tmp ,  dtype = float ) 

            fout.write('%15d%15.5f%15.5f\n' % (i * self.skip * self.timeStep , np.mean(tmp) , np.sqrt(np.var(tmp)/ len(tmp)   )   ))

        fout.close() 

if __name__ == '__main__':
    
    inputP = parseInput(sys.argv[1:])

    paraOpt = '-f -n -h -nsample  -skip -cutoff -v -o -ostat'.split(' ')

    helpdoc = ''\
              'Usage prog.py  -f trajNojump.xtc ; traj of nojump from g_traj  \n'\
              '               -n mol.ndx  ; relative diffusion between mol i and mol j , index start from 1 \n'\
              '               -nsample  10    ; number of random sample  \n'\
              '               -skip          ; number of frames to skip in survivial time  \n'\
              '               -cutoff          ;   \n'\
              '               -v          ; show progress if -v is set   \n'\
              '               -o out.xvg  ; output detail for each sampling; if not set do not print   \n'\
              '               -ostat out.stat  ; output statistics of nsamples \n'
    
    print_help(inputP, paraOpt , helpdoc )

    prob = Survival( inputP )

    prob.solv()
