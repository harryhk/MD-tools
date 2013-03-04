import sys , pdb
import numpy as np 
import re

#  input.sim : 
#  atomicdens: dens contribution from  each atom. First line Header, second line atomic contributions.  
#  index start:   
#  index end:  index =1 eq. column 1. Column 0 is z axis. Index end is included. 


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


inputP = parseInput(sys.argv[1:])
paraOpt = [ '-f', '-s', '-idx', '-z','-symm', '-h'  ]


for i in inputP:
    if i not in paraOpt or i=='-h':
        helpdoc ="Usage! ./prog.py\n"\
                "       -f input.sim  ; input sim file, first column: z axis, second-end : each type of atom  \n"\
                "       -s atomicdens ; First line atom type header ; second line dens for each atom type \n"\
                "       -idx index_s index_e  ; the sum from ith atom to jth atom. Closed brackets. Atom index start from 1\n"\
                "       -z z_min zmax ; set z range. if not set, use the whole data range. If -symm set, |z_min| = |z_max|\n"\
                "       -symm ; symmetrized around 0 \n"
        sys.stderr.write(helpdoc)
        sys.exit(1)

sim = np.array( [ i.split() for i in  open(inputP['-f']).readlines() if i[0]!='#' ] , dtype='float')
eN  = np.array( open(inputP['-s']).readlines()[1].strip().split()  , dtype='float')

sys.stderr.write( "Total number of atom types : %d \n" % len(eN) ) 

if len(eN) != sim.shape[1]-1:
    sys.stderr("sim file and dens file not match\n")
    sys.exit(1)

idx_s = int(inputP['-idx'][0])
idx_e = int(inputP['-idx'][1])

if idx_e > len(eN):
    sys.stderr("index_e out of range\n")
    sys.exit(1)

zmin = min( sim[:,0])
zmax = max( sim[:,0])
zflag = 0
if inputP.has_key('-z'):
    zflag =1 
    zmin = float(inputP['-z'][0])
    zmax = float(inputP['-z'][1])

symflag = 0
if inputP.has_key('-symm'):
    symflag = 1
    if zmin != -zmax:
        sys.stderr.write("symm set, zmin and zmax not differ by sign\n")
        sys.exit(1)

pdb.set_trace()
if zflag == 1 :
    sim = np.array([ i for i in sim if i[0]>=zmin and i[0] <= zmax ]  , dtype='float' )

i_b = 0 
i_e = sim.shape[0] - 1


if symflag == 1:
    while i_e - i_b >=1 :
        sim[i_b, 1:] , sim[i_e, 1:] = [ (sim[i_b, 1:] + sim[i_e, 1:])/2.0 ] * 2 

        i_b +=1 
        i_e -=1 
    

eN_tot = np.dot( sim[:, idx_s:idx_e+1] , eN[idx_s-1:idx_e])  # sim column start at 1. eN row start at 0

# write comment for the output 
print '# ' + ' '.join(sys.argv) 

for i in range(len(sim)):
	print "%15.6f%15.6f" %(sim[i,0] , eN_tot[i])


