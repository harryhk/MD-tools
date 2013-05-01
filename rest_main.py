# main file to modify topology for rest 
import sys
import rest_lipid_top, rest_force_field_top
import lnx_util

inputP = lnx_util.parseInput( sys.argv[1:] ) 
paraOpt = '-top  -bondedTop -nonbondedTop  -otop  -obonded  -ononbonded  -T1  -T2  -h'.split()

helpdoc = 'Usage! ./prog.py \n'\
          '       -top           top1  ; topology file for solutes that need rescale\n'\
          '       -bondedTop     top2  ; force field file that contains bonded and non_boned interactions\n'\
          '       -nonbondedTop  top3  ;\n'\
          '       -otop          outTop1 ; topology modified from top1\n'\
          '       -obonded       outTopbonded ; topology modified for bonded\n'\
          '       -ononbonded    outTopnonbonded ; topology modified for nonbonded\n'\
          '       -T1            323 k  ; \n'\
          '       -T2            600 k  ; gamma = beta_2 / beta_1      \n'\
          '                              ; charge: q_i = sqrt(gamma) *q_i ;\n'\
          '                              ; vdw   : epsilon_i = epsilon_i * gamma  ; i belongs to solute  ;\n'\
          '        -h                    ; print help doc \n'

lnx_util.print_help(inputP, paraOpt, helpdoc)

t1 = float( inputP['-T1'] )
t2 = float( inputP['-T2'] )
gamma = t1/ t2 

solute_top = rest_lipid_top.Topology( inputP['-top'], gamma)
bondedff = rest_force_field_top.BondedForceField(inputP['-bondedTop'], solute_top.bondMap, gamma )
nonbondedff = rest_force_field_top.NoneBondedForceField( inputP['-nonbondedTop'], solute_top.nameMap, gamma )


solute_top.display( open(inputP['-otop'], 'w' ) )
bondedff.display(open( inputP['-obonded'], 'w') )
nonbondedff.display(open( inputP['-ononbonded'], 'w') )


