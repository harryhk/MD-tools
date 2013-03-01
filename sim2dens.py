import sys , pdb
import numpy as np 

#  input.sim : 
#  atomicdens: dens contribution from  each atom. First line Header, second line atomic contributions.  
#  index start:   
#  index end:  index =1 eq. column 1. Column 0 is z axis. Index end is included. 


if len(sys.argv) != 5:
	print "Usage! ./prog.py input.sim atomicdens index_s index_e"
	exit(1)

sim = np.array( [ i.split() for i in  open(sys.argv[1]).readlines() if i[0]!='#' ] , dtype='float')
eN  = np.array( open(sys.argv[2]).readlines()[1].strip().split()  , dtype='float')

if len(eN) != sim.shape[1]-1:
    print "sim file and dens file not match"
    exit(1)

idx_s = int(sys.argv[3])
idx_e = int(sys.argv[4])

eN_tot = np.dot( sim[:, idx_s:idx_e+1] , eN[idx_s-1:idx_e])  # sim column start at 1. eN row start at 0

# write comment for the output 
print '# ' + ' '.join(sys.argv) 

for i in range(len(sim)):
	print "%10.6f%10.6f" %(sim[i,0] , eN_tot[i])


