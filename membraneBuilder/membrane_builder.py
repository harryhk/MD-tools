# this script used to build cholesterol and DOPC bilayer mixtures. 
# 

import numpy as np 
import sys , copy, random
import common.lnx_util as lu

def rotate_matrix(vec, theta): 
    ''' 
    return rotation matrix around axis vec about angle vec (in pi unit)
    '''
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    vec = vec / np.linalg.norm(vec)
    x, y, z = vec
    u_c = np.array( [ [0, -z, y],[z, 0, -x],[-y, x, 0] ] )
    return np.eye(3) * cos_theta + sin_theta * u_c + (1-cos_theta) * np.einsum('i,j', vec, vec)

def rotation_2vec(vec1, vec2): 
    '''
    return rotation matrix to rotate vec1 to vec2 
    '''
    vec3 = np.cross(vec1, vec2)
    vec3 = vec3 / np.linalg.norm(vec3)
    cos_theta = np.dot(vec1, vec2) / ( np.linalg.norm(vec1) * np.linalg.norm(vec2) )  
    sin_theta = np.sqrt( 1 - cos_theta * cos_theta )
    x,y,z = vec3
    u_c = np.array( [ [0, -z, y],[z, 0, -x],[-y, x, 0] ] )
    return np.eye(3) * cos_theta + sin_theta * u_c + (1-cos_theta) * np.einsum('i,j', vec3, vec3)


class Molecule(object):
    def __init__(self, lines): # molecules intitiated from lines in gro format
        self.numAtom = len(lines)
        self.resid = [  int(i[:5]) for i in lines ]
        self.resname =  [ i[5:10].strip() for i in lines  ]
        self.atomName = [ i[10:15].strip() for i in lines ]
        self.index = np.array( [ i[15:20] for i in lines] ,dtype=int )
        self.coor = np.array([ i[20:44].split() for i in lines ], dtype= float)

    def set_center(self):
        cog = np.mean(self.coor, 0)
        self.coor = self.coor - cog 

    def get_center(self):
        return np.mean( self.coor, 0 )

    def get_resname(self):
        return self.resname[0]
    
    def set_align(self, axis='z'):
        '''
        align the longest principal axis along z 
        '''
        inertia = np.dot(self.coor.transpose(), self.coor)
        inertia = np.eye(3)* inertia.trace() - inertia
        eigv , eigvec = np.linalg.eigh(inertia)
        if axis == 'z':
            rot_matrix = rotation_2vec( eigvec[:,0], (0,0,1)  )

        self.coor = np.dot( rot_matrix, self.coor.transpose() )
        self.coor = self.coor.transpose()
    
    def copy(self):
        '''
        return a deep copy of self object
        '''
        return copy.deepcopy(self)

    def set_rotate(self, axis, theta):
        rotM = rotate_matrix(axis, theta)
        self.coor = np.dot(rotM , self.coor.transpose())
        self.coor = self.coor.transpose()
        return self

    def set_trans(self, transvec):
        '''
        transfer the coor by transvec
        '''
        self.coor += transvec
        return self
      
    def set_resid(self, id):
        self.resid = np.ones_like(self.resid) * id % 100000 

    def set_index(self, idstart):
        self.index = ( self.index+ idstart - self.index[0] ) % 100000  # restart index if it is more than 5 digits

    def get_box(self):
        '''
        return the box dimension that can cover the molecule
        '''
        return tuple( self.coor.max(0) - self.coor.min(0) )

    def get_polar_box(self):
        '''
        return the box in polar coordinates
        '''
        r = np.max(  np.sqrt( self.coor[:,0] **2 + self.coor[:,1] **2  ) ) 
        z = self.coor[:,2].max() - self.coor[:,2].min()
        return (r, z)
    
    def get_coor_by_name(self, aName):
        '''
        return coordinates of a atom in molecule by name 
        '''
        return self.coor[ self.atomName.index(aName)  ]

    def print_gro(self, out):
        for resid, rN, aN, aIdx, coor in zip(self.resid, self.resname, self.atomName, self.index, self.coor ):
            print >> out, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" %( resid, rN, aN, aIdx, coor[0], coor[1], coor[2]) 


class System(object):
    def __init__(self):
        self.sys=[]
    
    
    def delete_filter(self, filter):
        '''
        delete molecule in sys if filter is true
        '''
        self.sys = [ i for i in self.sys if not filter(i) ]

        

    def extend(self, other):
        self.sys += copy.deepcopy(other.sys)

    def add_mol(self, mol):
        '''
        add a molecule into the sys
        '''
        assert isinstance(mol, Molecule)
        self.sys.append(mol.copy())

    def set_trans(self, transvec):
        for i in self.sys:
            i.set_trans(transvec)
        return self

    def get_center(self):
        c = 0 
        for i in self.sys:
            c += i.get_center()
        return c/ len(self.sys)
            
    
    def renumber(self):
        '''
        renumber the resid and index in each molecule
        '''
        index_start = 1
        for i in range(len(self.sys)):
            resid = i+1
            self.sys[i].set_resid(resid)
            self.sys[i].set_index(index_start)
            index_start += self.sys[i].numAtom
    
    def tot_atoms(self):
        return sum( [ i.numAtom  for i in self.sys ] )

    def tot_resids(self):
        return len(self.sys)

    def print_gro(self, out):
        print >> out , "system"
        print >> out , self.tot_atoms()
        self.renumber()
        for i in self.sys:
            i.print_gro(out)

class Monolayer(System):
    def __init__(self, nx, ny, dx, dy):
        System.__init__(self)
        self.nx , self.ny, self.dx, self.dy = ( nx, ny, dx, dy )
    
    def add_mols(self, mol, index, random_rotate=False):
        ''' 
        add mol into the index position 
        index is the list of 1d indexes in range(nx * ny)
        random_rotate if you want to rotate the adding mols
        '''
        index2grid =  lambda x: ( x % self.nx , x/self.nx )

        for i in index:
            idx_x , idx_y = index2grid(i)
            tmp_mol = mol.copy()
            # a random rotation around z axis 
            if random_rotate:
                tmp_mol.set_rotate( [0,0,1], 2* np.pi * random.random() )
            
            tmp_mol.set_trans( [idx_x * self.dx, idx_y * self.dy, 0] )
            System.add_mol(self, tmp_mol)
    


if __name__ == '__main__':

    inputP = lu.parseInput(sys.argv[1:])
    paraOpt = [ '-grid', '-chol', '-dopc', '-o', '-randomRotate'  ]

    
    helpdoc="Usage!  ./prog.py\n"\
            "        -grid 16 16     ; grid for monolayer\n"\
            "        -chol chol.base num ; only for monolayer\n"\
            "        -dopc dopc.base num\n"\
            "        -randomRotate   dopc chol    ; set it when you want to a random rotate on each molecule. You may need a much larger box\n"\
            "        -o    out.gro      \n"

    
    lu.print_help(inputP, paraOpt, helpdoc)

    # build dopc gro 
    dopc = Molecule( open(inputP['-dopc'][0]).readlines() )
    chol = Molecule( open(inputP['-chol'][0]).readlines() )
    
    dopc.set_center()
    dopc.set_align()

    chol.set_center()
    chol.set_align()
    
    nDopc = int(inputP['-dopc'][1])
    nChol = int(inputP['-chol'][1])

    # build grid 
    nx, ny = [ int(i) for i in inputP['-grid'] ]
    
    print "Constructing bilayers on a %d x %d grid" %( nx, ny )
    print "DOPC dimension : %8.3f%8.3f%8.3f" %( dopc.get_box() ) 
    print "DOPC dimension in polar coor: %8.3f%8.3f" %( dopc.get_polar_box() ) 
    print "CHOL dimension : %8.3f%8.3f%8.3f" %( chol.get_box() )
    print "CHOL dimension in polar coor: %8.3f%8.3f" %( chol.get_polar_box() )
    dx , dy = [  float(i) for i in  raw_input("Input Grid spacing dx, dy: ").split()  ]
        
    

    n_tot = nx * ny 
    # build upper layer 
    bilayer_u = Monolayer(nx, ny, dx, dy)
    chol_index = random.sample( range(n_tot), nChol )
    
    # rotate flag
    rand_rot_chol, rand_rot_dopc = (False, False) 
    try:
        rand_rot_chol = 'chol' in inputP['-randomRotate']
        rand_rot_dopc = 'dopc' in inputP['-randomRotate']
    except KeyError:
        pass
    
    bilayer_u.add_mols( dopc, list( set(range(n_tot)) - set(chol_index) ), rand_rot_dopc )
    bilayer_u.add_mols( chol, chol_index , rand_rot_chol)

    bilayer_l = Monolayer(nx, ny, dx, dy)
    chol_index = random.sample( range(n_tot), nChol )
    
    bilayer_l.add_mols( dopc.copy().set_rotate([1,0,0], np.pi), list( set(range(n_tot)) - set(chol_index) ), rand_rot_dopc )
    bilayer_l.add_mols( chol.copy().set_rotate([1,0,0], np.pi), chol_index , rand_rot_chol)

    bilayer_seq = float( raw_input("Separation between monolayer: (2.8) ")  )
    bilayer_u.set_trans([0,0, bilayer_seq])
    bilayer_u.extend(bilayer_l) 


    box_z = float( input("input box z (nm): ")  ) 
    bilayer_u.set_trans([0,0, box_z/2.0 - bilayer_u.get_center()[2]])

    fout = open(inputP['-o'], 'w')
    bilayer_u.print_gro(fout)
    print >> fout , '%5.1f%5.1f%5.1f' % ( nx* dx, ny*dy , box_z) 
