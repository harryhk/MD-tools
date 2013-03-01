import sys, lnx
import numpy as np 

if len(sys.argv)!= 3:
	sys.stdout.write("Usage: ./prog.py fin cut\n")
	sys.exit(1)

data= lnx.xvg(sys.argv[1])
cut = float( sys.argv[2] )

data = [ i for i in data if i[0] >= -cut and i[0] <= cut ]

for i in data:
	sys.stdout.write( '%15.4f' % i[0] + ''.join( [ '%15.6f' % ( j)   for j in i[1:] ]) + '\n'  )
 

