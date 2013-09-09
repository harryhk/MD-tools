import sys
from common.lnx_util import parseInput, print_help 

d={}
l=[]

inputP = parseInput(sys.argv[1:])
paraOpt = ['-f', '-h']

helpdoc="Usage!  ./prog.py\n"\
        "        -f  run.gro  ; input gro file. Should be gro file because the program strip the first 2 lines and last line for the box\n"\
        "        -h           ; print help document\n"

print_help(inputP, paraOpt, helpdoc)

gro = open(inputP['-f']).readlines()

for i in gro[2:-1]:
    atom_name = i[10:15].strip()
    atom_indx = int( i[15:20] ) 
    if atom_name not in d :
        d[atom_name] = [ atom_indx  ]
        l.append(atom_name)
    else:
        d[atom_name].append(atom_indx)

for i in l:
    print "[ %s ]" % i
    print "".join( [ ("%s " % j) for j in d[i] ] ) 

    
