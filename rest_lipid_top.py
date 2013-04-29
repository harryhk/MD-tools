# this file modify the lipid.itp file for rest 

import sys, re

debug = 1

def typeA2B(c):
    if not c.isalpha():
        return c

    d = c.upper()
    if c.islower():
        return chr( ord('A') + ( ord(d) - ord('A')  + 3) % 26 ).lower()
    else:
        return chr( ord('A') + ( ord(d) - ord('A')  + 3) % 26 )



class Atom(object):
    header = ( 'nr', 'type' , 'resnr' , 'residue' , 'atom' , 'cgnr' , 'charge' , 'mass', 'typeB', 'chargeB', 'massB' )
    header_string = '%8s' * len(header) % header
    format_string = "%8d%8s%8d%8s%8s%8d%10g%10g%8s%10g%10g"
    
    def __init__(self, line, scale = 1.0):
        '''
        use a single line to initialize Atom object
        line should not be the comment line in topology files
        '''
        line = re.sub( re.compile(';.*'), '', line  )  # line[:line.index(';')]
        self.nr , self.type , self.resnr , self.residue , self.atom , self.cgnr , self.charge , self.mass  =\
         [ _i(_j) for _i, _j in zip( [ int, str, int, str, str, int, float, float ] , line.strip().split() ) ]
        
        self.massB = self.mass
        self.typeB =''.join(map(typeA2B, self.type))
        self.chargeB = self.charge * scale

    def __repr__(self):
        return Atom.format_string % ( self.nr, self.type, self.resnr, self.residue, self.atom, self.cgnr, self.charge, self.mass, self.typeB, self.chargeB, self.massB)
    
    def display_typeAB(self):
        return Atom.format_string % ( self.nr, self.type, self.resnr, self.residue, self.atom, self.cgnr, self.charge, self.mass, self.typeB, self.chargeB, self.massB)


class Bonded(object):
    
    @staticmethod 
    def header_string(x):
        return  '%8s' * len(x) % x 
    
    def __init__(self, line, index):
        '''
        bonded interactions topology base class 
        index indicate the index for bonded parameters 
        items before index are int 
        '''
        line = re.sub( re.compile(';.*'), '', line  ).strip().split()   #  line[: line.index(';')].strip().split()
        self.data = map(int, line[:index] )
        self.par_A = line[index] if index == len(line)-1 else ''
        # parameter type should be 1 elements
        #assert len(self.par_A) <= 1
        #self.par_B = ''.join( map(typeA2B , self.par_A) )
        self.par_B = [ i for i in self.par_A ] 
        if self.par_B:
            self.par_B[0] = 'h'
            self.par_B = ''.join(self.par_B)

        
        self.index = index

    def get_string_AB(self):
        '''
        return string include type A and type B
        '''
        return ( '%8d' * self.index + '%12s' * 2)  % ( tuple(self.data) + (self.par_A, self.par_B) )
        

class Bond(Bonded):
    header = ('ai', 'aj', 'funct', 'par_A', 'par_B')
    header_string = Bonded.header_string(header) 
    
    def __init__(self, line):
        '''
        use a single line to initialize Bond object
        line should not be the comment line in topology files
        '''
        Bonded.__init__(self, line, 3)
        
    def __repr__(self):
        return self.get_string_AB() 

class Pair(Bonded):
    header = ('ai', 'aj', 'funct')
    header_string = Bonded.header_string(header) 
    
    format_string = '%8d' *3 

    def __init__(self, line):
        Bonded.__init__(self, line , 3)

    def __repr__(self):
        tmp = tuple(self.data) 
        return Pair.format_string % tmp 

class Angle(Bonded):
    header = ('ai', 'aj', 'ak', 'funct', 'par_A', 'par_B')
    header_string = Bonded.header_string(header)

    def __init__(self, line):
        Bonded.__init__(self, line , 4)

    def __repr__(self):
        return self.get_string_AB()

class Dihedral(Bonded):
    header = ('ai', 'aj', 'ak', 'al', 'funct', 'par_A', 'par_B' )
    header_string = Bonded.header_string(header)

    def __init__(self, line):
        Bonded.__init__(self, line, 5)

    def __repr__(self):
        return self.get_string_AB()

class Moleculetype(object):
    header = ( 'Name', 'nrexcl' )
    header_string = '%8s' *2 % header

    def __init__(self, line):
        self.line = re.sub( re.compile(';.*'), '', line  )  #line[: line.index(';')].strip().split()

    def __repr__(self):
        return self.line 

class MoleculeTop(object):
    flags = ( 'moleculetype', 'atoms', 'bonds', 'pairs', 'angles', 'dihedrals'  )
    top_class = { 'moleculetype' : Moleculetype,
                  'atoms'        : Atom,
                  'bonds'        : Bond,
                  'pairs'        : Pair,
                  'angles'       : Angle,
                  'dihedrals'    : Dihedral
                }

    def __init__(self, lines):
        '''
        lines to initialize a complete molecule 
        first line must start with [ moleculetype ]
        '''
        
        self.top={}
        self.flag = None 

        for line_num, line in enumerate(lines):
            if self.set_flag(line):
                continue
            else:
                try:
                    self.top[self.flag].append( MoleculeTop.top_class[ self.flag ](line)  ) 
                except KeyError:
                    self.top[self.flag]= [ MoleculeTop.top_class[ self.flag ](line)  ]

    def set_flag(self, line):
        if line[0] =='[':
            self.flag = line.strip('[] ')
            return True
        else:
            return False
    
    def atom_types(self):
        '''
        return a set of all atomtypes in the molecule
        '''
        return set( [ i.type for i in  self.top['atoms'] ] )
    
    def display(self, fout):
        for i in MoleculeTop.flags:
            # print header string 
            print >> fout, '[ ' +  i  + ' ]'
            print >> fout, '; '+ MoleculeTop.top_class[i].header_string

            # sub class
            for j in self.top[i]:
                print >> fout, repr(j)

            print >> fout 
        

class Topology(object):
    def __init__(self, filename):
        '''
        read the whole topology file at once 
        strip spaces at the begining and end of each line
        delete lines starts with ;
        split lines into complete moleculetype
        '''
        self.top = []

        lines = [ _i.strip() for _i in open(filename).readlines() ]
        lines = [ _i for _i in lines if len(_i) > 0 and  _i[0] != ';' ]
        mole_idx = [ i for i, j in enumerate(lines) if j.strip('[] ') == 'moleculetype' ]
        mole_idx.append(len(lines)) 

        for idx_start, idx_end  in zip( mole_idx[:-1], mole_idx[1:] ): 
            self.top.append( MoleculeTop(lines[idx_start: idx_end])  )

    def display(self, fout):
        for i in self.top:
            i.display(fout)

    def all_atom_Atypes(self):
        '''
        return a set of all atom Atypes 
        '''
        tmp = set()
        for i in self.top:
            tmp |= i.atom_types()
        return tmp 

            
            

    



if __name__ == '__main__':
    top = Topology(sys.argv[1])
    top.display(sys.stdout)
    print top.all_atom_Atypes()
