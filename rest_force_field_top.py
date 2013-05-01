# this file modify forcefield.itp (bonded and nonbonded) for rest

import re, math
import rest_lipid_top 

headerString = lambda x : '%8s' * len(x) % x 

class Atomtype(object):
    header = tuple('name     mass     charge     ptype     C6     C12'.split() )
    header_string = headerString(header) 
    format_string = '%5s' + '%10.3f' *2 + '%10s' + '%15.7E' * 2 

    def __init__(self, line):
        '''
        use a single line to initialize Atomtypes object
        line should not be the comment line 
        '''

        line = re.sub( re.compile(';.*'), '', line )
        self.name, self.mass, self.charge, self.ptype, self.c6, self.c12 = \
        [ _i(_j) for _i, _j in zip([str, float, float, str, float, float], line.strip().split() ) ]

        self.Bflag = False # if it has Btype 
    
    def __eq__(self, other):
        return self.name == other.name

    def typeB(self, nameMap, scale):
        '''
        nameMap store the oldname--> newname that requires to be changed 
        the names that need change are determined from rest_lipid_top.py
        !!! The scaling only support C6 , C12 format !!!
        '''
        assert isinstance(nameMap, dict)
        if self.name in nameMap:
            self.nameB = nameMap[self.name]
            self.c6B   = self.c6 * scale
            self.c12B  = self.c12 * scale
            self.Bflag = True

    def get_string(self):
        return Atomtype.format_string % (  self.name, self.mass, self.charge,  self.ptype, self.c6, self.c12  )

    def get_stringB(self):
        return Atomtype.format_string % (  self.nameB, self.mass, self.charge,  self.ptype, self.c6B, self.c12B  )


class NoneBonded(object):
    def __init__(self, line):
        line = re.sub( re.compile(';.*'), '', line )
        self.i, self.j , self.func, self.c6, self.c12 =\
        [ _i(_j) for _i, _j in zip([str, str, int, float, float], line.strip().split() ) ]
        self.Bflag = False

    def typeB(self, nameMap, scale):
        '''
        scale noneBonded parameters for typeB
        '''
        assert isinstance(nameMap, dict)
        i_inMap = self.i in nameMap
        j_inMap = self.j in nameMap

        if i_inMap or j_inMap:
            if i_inMap and (not j_inMap):
                self.iB = nameMap[self.i]
                self.jB = self.j
                self.c6B = self.c6 * scale
                self.c12B = self.c12 * scale 
             
            if (not i_inMap) and j_inMap: 
                self.iB = self.i
                self.jB = nameMap[self.j]
                self.c6B = self.c6 * scale
                self.c12B = self.c12 * scale 
            
            if i_inMap and j_inMap: 
                self.iB = nameMap[self.i]
                self.jB = nameMap[self.j]
                self.c6B = self.c6 * scale * scale
                self.c12B = self.c12 * scale * scale
            
            self.Bflag = True
            


    def get_string(self):
        '''
        default representing string
        '''
        return ('%8s' *2 + '%8d' + '%15.7E' *2 ) % (self.i, self.j, self.func, self.c6, self.c12 )

    def get_stringB(self):
        '''
        type B representing string
        '''
        return ('%8s' *2 + '%8d' + '%15.7E' *2 ) % (self.iB, self.jB, self.func, self.c6B, self.c12B )

class NoneBond(NoneBonded):
    header = tuple('i     j     func     C6      C12'.split())
    header_string = headerString(header)

    def __init__(self, line):
        NoneBonded.__init__(self, line)

    def __repr__(self):
        return self.get_string()

class Pairtype(NoneBonded):
    header = tuple('i     j     func        cs6        cs12'.split())
    header_string = headerString(header)

    def __init__(self, line):
        NoneBonded.__init__(self,line)

    def __repr__(self):
        return self.get_string()


class Bonded_type(object):
    def __init__(self, line):
        line = line.split()
        self.define , self.type  = ( line[0], line[1] ) 
        self.data  = map(float, line[2:] ) # self data has all bonded parameters 
        
        self.Bflag = False
        
    def _typeB(self, typeB, bondMap,  scale, indexes):
        '''
        bondMap has all the bondtype that need typeB
        index has the location(s) in self.data that needs scale, should be interable
        '''
        if self.type in bondMap:
            self.dataB = self.data[:]
            self.typeB = typeB 
            for j in indexes:
                self.dataB[j] = self.data[j] * scale
            self.Bflag = True
        
    def set_typeB(self, bondMap, scale):
        '''
        use bondMap to decide how to set scale
        internally call _typeB
        bondMap is a dict
        typeA : (typeB, funct, 'bonds or angles and etc')
        '''
        assert isinstance(bondMap, dict)
        if self.type in bondMap:
            typeB_name = bondMap[self.type][0]
            bond_type = ( bondMap[self.type][1], bondMap[self.type][2] )
            if bond_type == (2, 'bonds'):
                index = (1,)
            elif bond_type == (2, 'angles'):
                index = (1,)
            elif bond_type == (1, 'dihedrals'):
                index = (1,)
            elif bond_type == (2, 'dihedrals'):
                index = (1,)
            elif bond_type == (3, 'dihedrals'):
                index = (0,1,2,3,4,5)
            else:
                raise Exception("Unknown Bonded type")

            self._typeB(typeB_name, bondMap,  scale, index)


    def get_string(self):
        return  ('%-10s' * 2 + '%15g' * len(self.data) )  % ( (self.define,) +  (self.type,) + tuple(self.data) ) 

    def get_stringB(self):
        return ( '%-10s' * 2 + '%15g' * len(self.data) ) % ( (self.define,) +  (self.typeB,) + tuple(self.dataB) ) 





class BondedForceField(object):
    def __init__(self, filename, bondMap, gamma):
        '''
        Create bonded force field include typeB 
        change all lines start with #define and leave every other lines unchanged 
        '''
        lines = [ i.strip() for i in open(filename).readlines() ]
        self.top=[]
        for i in lines:
            if i[:7] == '#define':
                tmp = Bonded_type(i)
                tmp.set_typeB(bondMap, gamma)
                self.top.append(tmp.get_string() ) # add typeA 
                if tmp.Bflag:
                    self.top.append( tmp.get_stringB() ) # add typeB if existed
            else:
                self.top.append(i)
        
    def display(self, fout):
        for i in self.top:
            print >> fout , i


class NoneBondedForceField(object):
    def __init__(self, filename, nameMap, gamma):
        '''
        Create nonbonded force field include typeB
        Strip all comment and empty lines
        '''
        lines = [ i.strip() for i in open(filename).readlines() ]
        lines = [ i for i in lines if len(i) >0 and i[0] != ';' ]

        self.top={}
        self.type_class = { 'atomtypes' : Atomtype ,
                            'nonbond_params' : NoneBond ,
                            'pairtypes'  :  Pairtype 
                          }
        self.nameMap = nameMap
        self.gamma = gamma

        # find index for [ atomtypes ], [ nonbond_params ], [ pairtypes ]
        types_idx = [ i for i , j in enumerate(lines) if j[0]=='['  ]
        types_idx.append(len(lines) )

        for i in range(len(types_idx) -1 ):
            tmp_type = lines[types_idx[i] ].strip('[] ')
            assert tmp_type in self.type_class
            
            for line in lines[types_idx[i]+1 : types_idx[i+1]]:
                try:
                    self.top[tmp_type].append( self.type_class[tmp_type](line)  )
                except KeyError:
                    self.top[tmp_type] = [ self.type_class[tmp_type](line) ]
        
        # self consistence check that typeB names does not collide with original type names 
        assert len( set( nameMap.values()  ) & set( [j.name for j in self.top['atomtypes'] ] ) )== 0


    def display(self, fout):
        for i in self.top:
            print >> fout , '[ '+ i + ' ]'
            print >> fout, '; ' + self.type_class[i].header_string

            # subclass
            for j in self.top[i]:
                print >> fout, j.get_string()
                if i == 'atomtypes':
                    j.typeB( self.nameMap, self.gamma  )
                else:   
                    j.typeB( self.nameMap, math.sqrt(self.gamma) ) 
                
                if j.Bflag:
                    print >> fout , j.get_stringB()

        print >> fout 

                    



        
