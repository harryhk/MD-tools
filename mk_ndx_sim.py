import sys

d={}
l=[]
for i in sys.stdin.readlines():
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

    
