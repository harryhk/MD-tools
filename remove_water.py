import sys
import membrane_builder as mb
import lnx_util as lu

if __name__=='__main__':
    
    inputP = lu.parseInput(sys.argv[1:])
    paraOpt = ['-f', '-o', '-z']

    helpdoc = "Usage!  ./prog.py\n"\
              "   -z zmax  zmin   ; delete water between zmax and zmin\n"\
              "   -f input.gro    ; input membrane gro with water inside membrane\n"\
              "   -o out.out      ; output gro file with water removed\n"
    
    lu.print_help(inputP, paraOpt, helpdoc)

    # read gro file 
    fin = open( inputP['-f'] ).readlines()
    gro_start = 2
    gro_end   =  len(fin) 
    
    i_pre = gro_start

    sys = mb.System()


    for i in range(gro_start, gro_end):
        if fin[i][0:10].strip() == fin[i_pre][0:10].strip():
            i += 1 
        else:
            sys.add_mol( mb.Molecule( fin[i_pre: i] ) )
            i_pre = i 
            i += 1
    
    zmax, zmin = ( float(i) for i in inputP['-z'] )

    
    sys.delete_filter( lambda i: i.get_resname() == "SOL" and i.get_center()[2] < zmax and i.get_center()[2] > zmin  )
    

    fout = open(inputP['-o'], 'w')
    sys.print_gro(fout)
    print >> fout , fin[-1].strip('\n')
        
    
