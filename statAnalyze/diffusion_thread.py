# this code calculates the < \delta x(t) ^ 2 > between molecule i and molecule j 

import xdrfile as xdrfile 
import numpy as np 
import sys , copy , Queue, threading 
import multiprocessing as mp
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
        # self.msdAllPair : dictionary stores that msd of all pairs 

        # if -thread is set , we need to multi threading .
        # self.flagThread : if it is set, we do multi thread 
        # self.flagExitThread : exit flag for all thread when self.queue is empty  
        # self.queueLock  : lock for the Queue 
        # self.queue      : queue for ( i, j  ) mol pairs 
        # self.threads    : list of all threads and used for join 
        # self.nThreads   : number of threads 


        trajFn = inputP['-f']
        molFn  = inputP['-n']


        # read traj 
        tmpx = []
        # tmpbox = []
        tmpT1 = -1 
        tmpT2 = -1 
        
        self.msdAllPair = {} 

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
        
        # part for multi thread 
        if inputP.has_key('-thread'):
            self.flagThread = True 
        else:
            self.flagThread = False 

        
        if self.flagThread:
            
            # all the things to set thread environment
            self.queueLock = threading.Lock() 
            self.queue = Queue.Queue() 

            self.queueLock.acquire()
            for  i in self.molPair:
                self.queue.put( i ) 

            self.queueLock.release()
            self.threads = []

            self.nThreads = int( inputP['-thread'] )
            self.flagExitThread = False


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
        c=0
        tmpPer = 1
        nc = float( len(self.molPair) )
        for i in self.molPair:
            self.msdAllPair[i] = self.msd(i[0], i[1] ) 
            c+=1
            if c/ nc > tmpPer / 100.0:
                print "Progress: %d " % tmpPer + '%'
                tmpPer += 1 

        
        self.print_msd( )
    
    def print_msd(self):
        timeStep =  [ i * self.timeStep  for i in self.deltaF ] 

        # print header 
        headerstr = '%15s' %('# time') + ''.join( [ '%15s' %   str(i)+'-'+str(j)   for i, j in self.molPair  ] ) 
        print >> self.fout , headerstr
        
        for idx in range( len( timeStep)   ):
            self.fout.write('%15.5f' % timeStep[idx])

            for i in self.molPair:
                self.fout.write('%15.5f' % self.msdAllPair[i][idx] )

            self.fout.write('\n')


    def allPair_msd_thread(self):
        
        # create threads 
        for i in range( self.nThreads):
            thread = MSD_Thread( i, self)
            thread.start() 
            self.threads.append(thread)

        while not self.queue.empty():
            pass

        self.flagExitThread = True 

        for t in self.threads:
            t.join()

        #print "exit main thread "

        self.print_msd()

    def solv(self):
        
        if self.flagThread:
            self.allPair_msd_thread()
        else:
            self.allPair_msd()

class MSD_Thread(threading.Thread ):
    
    def __init__(self, id , prob  ):
        #  prob should be diffusion class  

        threading.Thread.__init__(self)
        self.threadId = id 
        self.prob = prob 

    def run(self) :
        # msdDict is a dictionary that stores all pairwise msd 
        while not self.prob.flagExitThread:
            
            self.prob.queueLock.acquire()
            if not self.prob.queue.empty():
                i = self.prob.queue.get() 
                self.prob.queueLock.release() 
                #print "Thread %d , pair %d %d" % (self.threadId, i[0], i[1])
                self.prob.msdAllPair[i] = self.prob.msd( i[0], i[1] )
            else:
                self.prob.queueLock.release()

        

if __name__ == '__main__':
    
    inputP = parseInput(sys.argv[1:])

    paraOpt = '-f -n -o -h -skip -thread'.split(' ')

    helpdoc = ''\
              'Usage prog.py  -f trajNojump.xtc ; traj of nojump from g_traj  \n'\
              '               -n mol.ndx  ; relative diffusion between mol i and mol j , index start from 1 \n'\
              '               -skip  1[defualt]    ; skip frames in delta T \n'\
              '               -thread  n           ; use n threads , beware that python with threads is pointless because of global interpreter lock \n'\
              '               -o out.xvg  ; \n'
    
    print_help(inputP, paraOpt , helpdoc )

    prob = Diffusion( inputP )

    prob.solv()
