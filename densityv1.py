# this code is designed to calculate the density of lipid bilayers specifically 
# it automatically assume that bilayer normal is along z axis, which is usually true. 
# bilayer center of mass is determined from a run.mass file which has following format:
#     atom index        mass 
#        1               15.035
#        2               15.035
#      ....              ...
# this file has all the atoms that in the bilayer 

# two ways can be used to calculate density. 
# 1. The absolutely position. zmin = 0 in this case and z value in the density histogram is the absoutely position in the box . (not default, set it with -absolute)
# 2. Use bilayer center as reference 0. This requires we to supply a reference atom index that belongs to the bilayer. We calcualte in the following ways:
#     first, center the system to the reference atom. 
#     second, define the center of mass of the bilayer as 0. 
# use -c atomindex to provide the reference atom. 

# all index starts from 1. 
# -f run.xtc 
# -s run.mass # basically the mass of all dopc atoms
# -n run.ndx 
# -o out.xvg 
# -b 200   # frame unit in xtc. if you write xtc every 500 steps and each step is 0.002 ps, then each frame is 1 ps. 
# -step  2  # calcualte every 2 frames in xtc file
# -e 400  (optional)
# -sl 200 
# -ng 1 (number of groups )
# -zrange -10 nm 10 nm  
# -outInAm 
# -h   # print help information 
# -c 508  ; index atom that used as center to remove pbc, if not set the absolute position is used.  

import sys, pdb 
import xdrfile as xdrfile 
import numpy as np
import re

def parseInput(l):
    inputP = {}
    
    flag = re.compile('-\D+')
    indx = [ i for (i, n ) in enumerate(l) if flag.match(n)   ]
    indx.append(len(l))
    
    for i in range( len(indx) -1 ):
        inputP[l[indx[i]]] = l[ indx[i]+1: indx[i+1] ]
        if len(inputP[l[indx[i]] ] ) ==1 :
            inputP[ l[indx[i]]] = inputP[l[indx[i]]][0]

    return inputP

def parseNdx(inputP):
    fin = [ i.strip() for i in   open(inputP['-n']) ]
    fin = [ i for i in fin if len(i) > 0  ]

    indx = [ i for (i,n) in enumerate(fin) if n[0] == '['  ]

    if len(indx) != int( inputP['-ng'] ) :
        sys.stderr.write('Index file problem')
        sys.exit(1)
    
    index_l = []
    for i in range(len(indx) ):
        if i != len(indx)-1:
            l = ' '.join( fin[ indx[i]+1 : indx[i+1] ] )
            index_l.append( l.split()  ) 
        else:
            l = ' '.join( fin[ indx[i]+1 :  ] )
            index_l.append(  l.split()  )
    
    for i in range( len(indx)  ):
        print fin[ indx[i]] , len( index_l[i] ) 
    
    return index_l 


class Mol(object):
    def __init__(self, index, mass= None, nbins = -1 ):
        self.index = np.array( index, dtype = int  )
        self.index -= 1 
        self.mass = np.array( mass, dtype =float  )

        if nbins > 0 :
            self.bins = np.zeros( nbins )

    def get_frame(self, x):  # only with z dimension 
        self.pos = np.array( [  i[2]  for i in  x[self.index] ] )
    
    def get_com(self):
        return  sum( self.pos * self.mass ) / sum(self.mass)
    
    def update_dens(self,  zs, binsize, box):
        (bx,by,bz) = np.diagonal(box)
        binIdx  =   np.floor( ( self.pos  - zs ) / binsize ).astype(int)
        dens = 1.0/(binsize * bx * by )
        for i in binIdx:
            self.bins[ i ] += dens
    
    def update_sym_dens(self, zc, zs, binsize, box   ):
        (bx, by, bz ) = np.diagonal(box) 
        dist = self.pos - zc
        dist = dist - np.round( dist / bz  ) * bz 
        binIdx = np.floor( ( dist - zs)/ binsize).astype(int)
        dens = 1.0/(binsize * bx * by )
        for i in binIdx:
            self.bins[ i ] += dens
        

def shift_z( oldx , cIdx , boxz   ): # shift system according to center atom cIdx; only center along z 
    cIdx -=1     # index start from 0 
    zc = oldx[cIdx][2]
    oldx[:,2] -= zc
    oldx[:,2] =  oldx[:,2] - np.round( ( oldx[:,2]  ) / boxz ) * boxz 
    


if __name__ == '__main__':
    
    inputP = parseInput( sys.argv[1:] )

    # check input 
    paraOpt = [ '-f', '-s', '-n', '-o', '-b', '-e', '-sl', '-ng', '-zrange', '-outInAm',  '-h' , '-c', '-step' ]

    for i in inputP:
        if i not in paraOpt:
            sys.stderr.write('Input parameter no recognized. Check with -h option \n')
            sys.exit(1)


    if inputP.has_key('-h'):
        helpdoc = '' \
                  'Usage!  ./prog.py   -f run1.xtc run2.xtc ...  \n' \
                  '                    -s bilayer.mass ; mass file of all bilayer atoms   \n' \
                  '                    -n run.ndx  ; gromacs index file of all the groups to compute density  \n' \
                  '                    -o out.xvg   \n' \
                  '                    -b 200      ; beginning of frame (ps) \n' \
                  '                    -e 400      ; ending of frame (ps)    \n' \
                  '                    -step  2    ; calcualte every 2 frames xtc file\n'\
                  '                    -sl 200     ; number of slices in z   \n' \
                  '                    -ng 1       ; number of groups to compute density. Should be consistent with run.ndx \n' \
                  '                    -zrange  -10 (nm) 10 (nm) ; range of z \n' \
                  '                    -outInAm    ; output z in Am unit instead of nm \n' \
                  '                    -c 508      ; index atom that used as center to remove pbc.\n'\
                  '                    -h          ; print help information\n'
                

        sys.stdout.write(helpdoc)
        sys.exit(1)



    # setup histogram 
    nBins = int( inputP['-sl'] ) 
    zRange_min = float( inputP['-zrange'][0]  )
    zRange_max = float( inputP['-zrange'][1]  )
    binSize = (zRange_max - zRange_min) / nBins;

    histogramBinCenter = np.linspace( zRange_min, zRange_max, nBins+1)
    histogramBinCenter = np.array( [ (i+j)/2.0 for i,j in zip( histogramBinCenter[0:-1] , histogramBinCenter[1:]  )] , dtype =float )
    
    index = parseNdx(inputP)

    mol = []
    for i in index:
        mol.append( Mol(i, nbins = nBins )  )


    # create bilayer with mass
    (dopcIdx, dopcMass) =  zip(*[   i.strip().split()   for i in open( inputP['-s'], 'r') ] )
    dopc = Mol( dopcIdx, mass= dopcMass )
    
    # set time for calculation
    time_start = -1 
    time_end   = np.inf  # ps 
    time_step = 1
    if inputP.has_key('-b'):
        time_start = float( inputP['-b']  )
    if inputP.has_key('-e'):
        time_end = float( inputP['-e'] )
    if inputP.has_key('-step'):
        time_step = int ( inputP['-step'] )
   
    # set symmetry flag
    symmFlag = False
    if inputP.has_key('-c'):
        symmFlag = True
        cIdx = int( inputP['-c'] )

    frame_c =0
    global_time = 0 
    
    if type(inputP['-f']) is not list:
        inputP['-f'] = [ inputP['-f']  ]

    for xtc_f in inputP['-f']:
        

        for f in xdrfile.xdrfile(xtc_f) : 
            
            global_time += 1
            # print every 100 ps 
            if global_time % 100 == 0 :
                sys.stdout.write("reading frame : %10.3f frames\r" % (global_time) ) 
                sys.stdout.flush()

            if time_step  > 1:
                if global_time % time_step != 0 : continue


            if global_time > time_end:
                break 

            if global_time >= time_start :
                # shift system 
                if symmFlag:
                    shift_z(f.x, cIdx, f.box[2][2] )

                    # calculate bilayer center 
                    dopc.get_frame(f.x)
                    dopc_c = dopc.get_com()

                    for grp in mol:
                        grp.get_frame(f.x)
                        grp.update_sym_dens( dopc_c, zRange_min, binSize, f.box )
                else:
                    for grp in mol:
                        grp.get_frame(f.x)
                        grp.update_dens(zRange_min, binSize, f.box)

                
                frame_c += 1

    
    fout = open(inputP['-o'], 'w')
    if inputP.has_key('-outInAm'):
        histogramBinCenter *= 10.0
    
    pFmt = lambda x : '%10.4f' % x 
    pFmtAm = lambda x : '%10.6f' % x 
    
    # write comment for the output:
    fout.write( '# ' +' '.join(sys.argv) + '\n'  )

    for i in range( nBins ):
        l = pFmt( histogramBinCenter[i]  )
        for j in mol:
            if inputP.has_key('-outInAm'):
                l += pFmtAm( j.bins[i] /frame_c / 1000.0 )
            else:
                l += pFmt( j.bins[i] /frame_c )
        fout.write( '%s\n' % l  )



