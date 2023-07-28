import numpy as np
import sys
import argparse
import os

parser = argparse.ArgumentParser(description='')
parser.add_argument('-p', '--path', help='path to simulation folder - REQUIRED', required=True, type=str)
parser.add_argument('-ron', '--ron', help='imposed growth rate [monomers/second]', required=False, type=float, default=8.0)
parser.add_argument('-thyd', '--thyd', help='hydrolysis time [seconds]', required=False, type=float, default=5.0)
parser.add_argument('-rnuc', '--rnuc', help='imposed nucleation rate [filaments/second]', required=False, type=float, default=1.0)
parser.add_argument('-Kbond', '--Kbond', help='bond constant [kT/sigma2]', required=False, type=float, default=1000.0)
parser.add_argument('-Kbend', '--Kbend', help='bending constant [kT/sigma2]', required=False, type=float, default=800.0)
parser.add_argument('-Kobst', '--Kobst', help='bending constant of the obstacles [kT/sigma2]', required=False, type=float, default=100.0)
parser.add_argument('-tstep', '--tstep', help='simulation timestep [seconds]', required=False, type=float, default=0.001)
parser.add_argument('-runtime', '--runtime', help='simulation run time [seconds]', required=False, type=float, default=600.0)
parser.add_argument('-frate', '--frate', help='frame rate [seconds]', required=False, type=float, default=1.0)
parser.add_argument('-sd','--seed', help='random number generator seed', required=False, type=int, default=1234)
parser.add_argument('-L','--L', help='box size [sigma]', required=False, type=float, default=200)
parser.add_argument('-ICNfils','--ICNfils', help='number of filaments in Initial Conditions [0 or 1] -- if more than 1 then creates circular obstacle of that size; if less than 0 then creates filament across X = 0', required=False, type=int, default=0)
parser.add_argument('-arrest','--arrest', help='arrest treadmilling dynamics: turn shrinking off', required=False, action = 'store_true')
parser.add_argument('-fixlength','--fixlength', help='fix the length of the filaments', required=False, type=float, default=0)
parser.add_argument('-Xbias','--Xbias', help='only nucleate filaments along the positive X direction', required=False, action = 'store_true')
parser.add_argument('-attraction','--attraction', help='turn on interfilament attraction', required=False, action = 'store_true')
parser.add_argument('-saturate','--saturate', help='set a maximum number of particles in the system', required=False, type=int, default=0)

args = parser.parse_args()
gpath = args.path
ron = float(args.ron)
thyd = float(args.thyd)
rnuc = float(args.rnuc)
Kbond = float(args.Kbond)
Kbend = float(args.Kbend)
Kobst = float(args.Kobst)
tstep = float(args.tstep)
runtime = float(args.runtime)
frate = float(args.frate)
seed = int(args.seed)
L = float(args.L)
ICNfils = int(args.ICNfils)
arrest = args.arrest
fixL = float(args.fixlength)
Xbias = args.Xbias
attraction = args.attraction
saturate = int(args.saturate)

if saturate == 0:
    saturate = L*L+10

# Initialise numpy's RNG
np.random.seed(seed)

poff = 1.0
if arrest:
    poff = 0.0

r = os.system('mkdir %s'%(gpath))
r = os.system('mkdir %s/Reactions'%(gpath))
r = os.system('cp -r Reactions/*  %s/Reactions'%(gpath))

f = open('%s/info.txt'%(gpath), 'w')
f.write("path:\t\t%s\n"%(gpath))
f.write("ron [mons/s]:\t\t%.1f\n"%(ron))
f.write("thyd [s]:\t\t%.1f\n"%(thyd))
f.write("rnuc [fils/s]:\t\t%.1f\n"%(rnuc))
f.write("Kbond [kT/sigma2]:\t%.1f\n"%(Kbond))
f.write("Kbend [kT/sigma2]:\t%.1f\n"%(Kbend))
f.write("tstep [s]:\t\t%.5f\n"%(tstep))
f.write("runtime [s]:\t\t%.1f\n"%(runtime))
f.write("frate [s]:\t\t%.1f\n"%(frate))
f.write("seed:\t\t\t%d\n"%(seed))
f.write("box size:\t\t%.1f\n"%(L))
f.write("poff:\t\t%.1f\n"%(poff))
f.close()

if ICNfils == 1:
    alpha = 2*np.pi*np.random.random()
    TX = 0.5*np.cos(alpha+np.pi)
    TY = 0.5*np.sin(alpha+np.pi)
    HX = 0.5*np.cos(alpha)
    HY = 0.5*np.sin(alpha)
    f = open('%s/configuration.txt'%(gpath), 'w')
    f.write('''First line of this test
4 atoms
2 bonds
0 angles
4 atom types
1 bond types
2 angle types
''')
    f.write('%.1f %.1f xlo xhi\n'%(-L/2,L/2))
    f.write('%.1f %.1f ylo yhi\n'%(-L/2,L/2))
    f.write('''-4.25 0.25 zlo zhi

Masses

1 1
2 1
3 1
4 1


Atoms

1 1 2 %f %f 0.0
2 1 3 %f %f 0.0
3 0 4 0.0 -0.5 -2.0
4 0 4 0.0 0.5 -2.0


Bonds

1 1 1 2
2 1 3 4
'''%(TX, TY, HX, HY))
    f.close()
elif ICNfils == 0:
    f = open('%s/configuration.txt'%(gpath), 'w')
    f.write('''First line of this test
2 atoms
1 bonds
0 angles
4 atom types
1 bond types
2 angle types
''')
    f.write('%.1f %.1f xlo xhi\n'%(-L/2,L/2))
    f.write('%.1f %.1f ylo yhi\n'%(-L/2,L/2))
    f.write('''-4.25 0.25 zlo zhi

Masses

1 1
2 1
3 1
4 1


Atoms

1 0 4 0.0 -0.5 -2.0
2 0 4 0.0 0.5 -2.0


Bonds

1 1 1 2
''')
    f.close()
elif ICNfils > 1: ## RING POLYMER! - L = N = ICNfils
    N = ICNfils
    r = N/(2*np.pi)
    a = 1/r
    f = open('%s/configuration.txt'%(gpath), 'w')
    f.write('''First line of this test
%d atoms
%d bonds
%d angles
4 atom types
1 bond types
2 angle types
'''%(N+2,N+1,N))
    f.write('%.1f %.1f xlo xhi\n'%(-L/2,L/2))
    f.write('%.1f %.1f ylo yhi\n'%(-L/2,L/2))
    f.write('''-4.25 0.25 zlo zhi

Masses

1 1
2 1
3 1
4 1


Atoms

''')
    f.write("1 0 4 0.0 -0.5 -2.0\n")
    f.write("2 0 4 0.0 0.5 -2.0\n")
    for i in range(N):
        X = r*np.cos(i*a)
        Y = r*np.sin(i*a)
        f.write("%d 1 1 %f %f 0.0\n"%(3+i,X,Y))
    f.write('''

Bonds

1 1 1 2
''')
    for i in range(N-1):
        f.write("%d 1 %d %d\n"%(2+i,3+i,4+i))
    f.write("%d 1 %d %d\n"%(3+i,4+i,3))
    f.write('''

Angles

''')
    for i in range(N-2):
        f.write("%d 2 %d %d %d\n"%(i+1,3+i,4+i,5+i))
    f.write("%d 2 %d %d %d\n"%(i+2,4+i,5+i,3))
    f.write("%d 2 %d %d %d\n"%(i+3,5+i,3,4))
    f.close()
elif ICNfils < 0: ## LINE POLYMER! - cross the box at X = 0 w/ L = -ICNfls
    N = int(-ICNfils)
    natoms = N+2
    nbonds = N
    nangls = N-2
    if N == int(L):
        nbonds = N+1
        nangls = N
    f = open('%s/configuration.txt'%(gpath), 'w')
    f.write('''First line of this test
%d atoms
%d bonds
%d angles
4 atom types
1 bond types
2 angle types
'''%(natoms,nbonds,nangls))
    f.write('%.1f %.1f xlo xhi\n'%(-L/2,L/2))
    f.write('%.1f %.1f ylo yhi\n'%(-L/2,L/2))
    f.write('''-4.25 0.25 zlo zhi

Masses

1 1
2 1
3 1
4 1


Atoms

''')
    f.write("1 0 4 0.0 -0.5 -2.0\n")
    f.write("2 0 4 0.0 0.5 -2.0\n")
    Y0 = -N/2.0+0.5
    for i in range(N):
        f.write("%d 1 1 %f %f 0.0\n"%(3+i,0.0,Y0+i))
    f.write('''

Bonds

1 1 1 2
''')
    for i in range(N-1):
        f.write("%d 1 %d %d\n"%(2+i,3+i,4+i))
    if N == int(L):
        f.write("%d 1 %d %d\n"%(3+i,4+i,3))
    f.write('''

Angles

''')
    for i in range(N-2):
        f.write("%d 2 %d %d %d\n"%(i+1,3+i,4+i,5+i))
    if N == int(L):
        f.write("%d 2 %d %d %d\n"%(i+2,4+i,5+i,3))
        f.write("%d 2 %d %d %d\n"%(i+3,5+i,3,4))
    f.close()

f = open('%s/in.local'%(gpath), 'w')
f.write('''log                 log.txt
units               lj
dimension           2
atom_style          molecular
read_data           configuration.txt extra/bond/per/atom 5  extra/special/per/atom 20  extra/angle/per/atom 3
''')
f.write("variable            ron equal %.1f                                 # growth rate [monomers/s]\n"%(ron))
f.write("variable            tauhyd equal %.1f                              # hydrolysis time [seconds]\n"%(thyd))
f.write("variable            rnuc equal %.1f                                # nucleation rate [filaments/s]\n"%(rnuc))
f.write("variable            poff equal %.1f                                # shrinking reaction initiation probability (0 or 1)\n"%(poff))
f.write("variable            Kbond equal %.1f                               # bond constant [kT/sigma2]\n"%(Kbond))
f.write("variable            Kbend equal %.1f                               # bend constant [kT/sigma2]\n"%(Kbend))
f.write("variable            Kobst equal %.1f                               # bend constant [kT/sigma2]\n"%(Kobst))
f.write("variable            tstep equal %f                                 # simulation timestep size [seconds]\n"%(tstep))
f.write("variable            run_time equal %.1f                            # simulation run time [seconds]\n"%(runtime))
f.write("variable            frame_rate equal %.1f                          # dumping interval [seconds]\n"%(frate))
f.write("variable            seed equal %d                                  # random number generator seed\n"%(seed))
f.write("variable            maxatoms equal %d                              # maximum number of atoms allowed in the system\n"%(saturate))
f.write('''variable            condatoms equal "atoms >= v_maxatoms+2"
variable            pon equal (1-v_condatoms)                     # growing reaction initiation probability (0 or 1)

variable            kon atom ${ron}/10.0                          # growth probability
variable            knuc atom ${rnuc}/10.0                        # nucleation probability
variable            thyd equal ${tauhyd}/${tstep}                 # hydrolysis time [simulation steps]
variable            rstep equal 0.1/${tstep}                      # reaction interval [simulation steps]
variable            run_steps equal ${run_time}/${tstep}          # simulation run time [simulation steps]
variable            dump_time equal ${frame_rate}/${tstep}        # dumping interval [simulation steps]

variable            stab_steps equal 1

group               ghosts type 4

special_bonds       lj 1.0 1.0 1.0
bond_style          harmonic
bond_coeff          1 ${Kbond} 1.0

angle_style         harmonic
angle_coeff         1 ${Kbend} 180.0
angle_coeff         2 ${Kobst} 180.0

pair_style          hybrid/overlay zero 1.50 cosine/squared 1.50
pair_coeff          * * cosine/squared 0.00 1.00 1.10
pair_coeff          * * zero 1.50
''')
if attraction:
    f.write("pair_coeff          1 1 cosine/squared 10.00 1.00 1.50 wca\n")
else:
    f.write("pair_coeff          1 1 cosine/squared 1.00 1.00 1.00 wca\n")
f.write('''pair_coeff          1 2 cosine/squared 1.00 1.00 1.00 wca
pair_coeff          1 3 cosine/squared 1.00 1.00 1.00 wca
pair_coeff          2 2 cosine/squared 1.00 1.00 1.00 wca
pair_coeff          2 3 cosine/squared 1.00 1.00 1.00 wca
pair_coeff          3 3 cosine/squared 1.00 1.00 1.00 wca

''')
if attraction:
    f.write("neigh_modify        exclude molecule/intra all\n\n")
f.write('''# Nucleation reactions molecular templates
molecule            mPreNucleation Reactions/pre_Nucleation.txt
molecule            mPostNucleation Reactions/post_Nucleation.txt

# Growth reactions molecular templates
molecule            mPreDimerOn Reactions/pre_DimerOn.txt
molecule            mPostDimerOn Reactions/post_DimerOn.txt
molecule            mPreTrimerOn Reactions/pre_TrimerOn.txt
molecule            mPostTrimerOn Reactions/post_TrimerOn.txt
molecule            mPreOligomerOn Reactions/pre_OligomerOn.txt
molecule            mPostOligomerOn Reactions/post_OligomerOn.txt

# Shrink reactions molecular templates
molecule            mPreDimerOff Reactions/pre_DimerOff.txt
molecule            mPostDimerOff Reactions/post_DimerOff.txt
molecule            mPreTrimerOff Reactions/pre_TrimerOff.txt
molecule            mPostTrimerOff Reactions/post_TrimerOff.txt
molecule            mPreQuartomerOff Reactions/pre_QuartomerOff.txt
molecule            mPostQuartomerOff Reactions/post_QuartomerOff.txt
molecule            mPreOligomerOff Reactions/pre_OligomerOff.txt
molecule            mPostOligomerOff Reactions/post_OligomerOff.txt

variable            vMaskHBefore atom "type==2 || type==3"
group               HeadBefore dynamic all var vMaskHBefore every 1
variable            vMaskCBefore atom "type==1"
group               CoreBefore dynamic all var vMaskCBefore every 1
variable            vMaskSAZ atom "v_vSA == 0"
fix                 fSAZ all store/state 10 v_vMaskSAZ

fix                 freact all bond/react  stabilization yes AllAtoms 0.1  reset_mol_ids no         &
                ''')
if Xbias:
    f.write("react Nucleation all ${rstep} 0.900000 1.100000 mPreNucleation mPostNucleation Reactions/map_Nucleation.txt prob v_pon ${seed} stabilize_steps ${stab_steps} modify_create overlap 0.9 modify_create nuc xor         &\n")
else:
    f.write("react Nucleation all ${rstep} 0.900000 1.100000 mPreNucleation mPostNucleation Reactions/map_Nucleation.txt prob v_pon ${seed} stabilize_steps ${stab_steps} modify_create overlap 0.9 modify_create nuc yes         &\n")
f.write('''                react DimerOn all ${rstep} 0.900000 1.100000 mPreDimerOn mPostDimerOn Reactions/map_DimerOn.txt prob v_pon ${seed} stabilize_steps ${stab_steps} modify_create fit 1 modify_create overlap 0.9         &
                react TrimerOn all ${rstep} 0.900000 1.100000 mPreTrimerOn mPostTrimerOn Reactions/map_TrimerOn.txt prob v_pon ${seed} stabilize_steps ${stab_steps} modify_create fit 1 modify_create overlap 0.9         &
                react OligomerOn all ${rstep} 0.900000 1.100000 mPreOligomerOn mPostOligomerOn Reactions/map_OligomerOn.txt prob v_pon ${seed} stabilize_steps ${stab_steps} modify_create fit 1 modify_create overlap 0.9         &
                react OligomerOff all ${rstep} 0.900000 1.100000 mPreOligomerOff mPostOligomerOff Reactions/map_OligomerOff.txt prob ${poff} ${seed} stabilize_steps ${stab_steps}         &
                react QuartomerOff all ${rstep} 0.900000 1.100000 mPreQuartomerOff mPostQuartomerOff Reactions/map_QuartomerOff.txt prob ${poff} ${seed} stabilize_steps ${stab_steps}         &
                react TrimerOff all ${rstep} 0.900000 1.100000 mPreTrimerOff mPostTrimerOff Reactions/map_TrimerOff.txt prob ${poff} ${seed} stabilize_steps ${stab_steps}         &
                react DimerOff all ${rstep} 0.900000 1.100000 mPreDimerOff mPostDimerOff Reactions/map_DimerOff.txt prob ${poff} ${seed} stabilize_steps ${stab_steps}

variable            vMaskHeadType atom "type==3"
group               HeadMons dynamic all var vMaskHeadType every 1
variable            vMaskTailType atom "type==2"
group               TailMons dynamic all var vMaskTailType every 1
variable            vMaskHTType atom "type==2 || type==3"
group               HTMons dynamic all var vMaskHTType every 1
variable            vMaskAlive atom "type==1 || type==2 || type==3"
group               alive dynamic all var vMaskAlive every 1

variable            vR atom "gmask(bond_react_MASTER_group)"
variable            vHR atom "gmask(HTMons) * gmask(bond_react_MASTER_group)"
fix                 fHR all store/state 10 v_vHR
variable            vMaskFHR atom "f_fHR == 1"
group               gFHR dynamic all var vMaskFHR every 1
variable            vMaskVHR atom "v_vHR == 1"
group               gVHR dynamic all var vMaskVHR every 1
variable            vHRBefore atom "gmask(HeadBefore) * gmask(bond_react_MASTER_group)"
variable            vNHRB atom "(1-gmask(HeadBefore)) * gmask(bond_react_MASTER_group) * (1-gmask(CoreBefore))"
variable            vSAZT atom "(1-gmask(TailMons)) * f_fSAZ"
variable            v1T atom "v_vNHRB * v_vHR"
variable            vSI atom "f_fSI + v_v1T * (step - f_fSI)"
fix                 fSI all store/state 1 v_vSI
variable            vSA atom (step-f_fSI)
variable            vTailsTime atom v_vSA/v_thyd
variable            vTailsE atom exp(-v_vTailsTime)
variable            vTailsP atom 1.0-exp(-v_vTailsTime)

fix                 fLang all langevin 1.0 1.0 1.0 ${seed}
fix                 fNVE AllAtoms_REACT nve

variable            realtime equal step*${tstep}

dump                1 all custom ${dump_time} output.xyz id mol type x y z fx fy fz f_fSI v_vTailsTime
dump_modify         1 format line "%d %d %d %.2f %.2f %.2f %.2f %.2f %.2f %.1f %.1f"

thermo              ${dump_time}
compute_modify      thermo_temp dynamic/dof yes
thermo_style        custom step v_realtime temp pe ke etotal epair ebond eangle press vol density atoms

compute             cTempAll all temp
compute_modify      cTempAll dynamic/dof yes
compute             cTempAlive alive temp
compute_modify      cTempAlive dynamic/dof yes
compute             peratom all pe/atom
compute             peAll all reduce ave c_peratom
compute             peAlive alive reduce ave c_peratom
variable            vNumAlive equal count(alive)
variable            vNumHeads equal count(HeadMons)
variable            vNumTails equal count(TailMons)

compute             cBonds all property/local batom1 batom2 btype
compute             cBondDxys all bond/local engpot force dist

dump                2 all local ${dump_time} bonds.dump c_cBonds[*] c_cBondDxys[*]
dump_modify         2 format line "%f %f %f %.2f %.2f %.2f"

fix                 print_thermo alive ave/time ${dump_time} 1 ${dump_time} c_cTempAll c_cTempAlive v_vNumAlive v_vNumHeads v_vNumTails c_peAll c_peAlive file thermo.txt

fix                 twodim all enforce2d

timestep            ${tstep}
run                 ${run_steps}
''')
f.close()
