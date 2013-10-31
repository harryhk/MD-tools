import sys
from common.lnx_util import parseInput, print_help 

# bug : DOPC and TYR could share the same atom name N but we want to separate them into different index file 
# fix : for index name, we group atom name and residue name together  [ N DOPC ].

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
    res_name = i[5:10].strip()
    if res_name =='DOPC' or res_name == 'CL' or res_name == 'SOL':
        type = ( atom_name, res_name )
    else:  # speical case such one tat have several arginines. Different then by their index 
        res_name1 = i[0:10].strip()
        type = ( atom_name, res_name1 )

    atom_indx = int( i[15:20] ) 
    if type  not in d :
        d[ type ] = [ atom_indx  ]
        l.append( type )  
    else:
        d[ type ].append(atom_indx)

for i in l:
    print "[ %s %s ]" % ( i[0] , i[1] )
    print "".join( [ ("%s " % j) for j in d[i] ] ) 

    
