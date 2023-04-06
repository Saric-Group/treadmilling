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
parser.add_argument('-tstep', '--tstep', help='simulation timestep [seconds]', required=False, type=float, default=0.001)
parser.add_argument('-runtime', '--runtime', help='simulation run time [seconds]', required=False, type=float, default=600.0)
parser.add_argument('-frate', '--frate', help='frame rate [seconds]', required=False, type=float, default=1.0)
parser.add_argument('-sd','--seed', help='random number generator seed', required=False, type=int, default=1234)

args = parser.parse_args()
gpath = args.path
ron = float(args.ron)
thyd = float(args.thyd)
rnuc = float(args.rnuc)
Kbond = float(args.Kbond)
Kbend = float(args.Kbend)
tstep = float(args.tstep)
runtime = float(args.runtime)
frate = float(args.frate)
seed = int(args.seed)

r = os.system('mkdir %s'%(gpath))
r = os.system('mkdir %s/Reactions'%(gpath))
r = os.system('cp -r Reactions/*  %s/Reactions'%(gpath))
r = os.system('cp configuration.txt %s/configuration.txt'%(gpath))

f = open('%s/info.txt'%(gpath), 'w')
f.write("path:\t\t%s\n"%(gpath))
f.write("ron [mons/s]:\t\t%.1f\n"%(ron))
f.write("thyd [s]:\t\t%.1f\n"%(thyd))
f.write("rnuc [fils/s]:\t\t%.1f\n"%(rnuc))
f.write("Kbond [kT/sigma2]:\t\t%.1f\n"%(Kbond))
f.write("Kbend [kT/sigma2]:\t\t%.1f\n"%(Kbend))
f.write("tstep [s]:\t\t%.1f\n"%(tstep))
f.write("runtime [s]:\t\t%.1f\n"%(runtime))
f.write("frate [s]:\t\t%.1f\n"%(frate))
f.write("seed:\t\t%d\n"%(seed))
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
f.write("variable            Kbond equal %.1f                               # bond constant [kT/sigma2]\n"%(Kbond))
f.write("variable            Kbend equal %.1f                               # bend constant [kT/sigma2]\n"%(Kbend))
f.write("variable            tstep equal %f                                 # simulation timestep size [seconds]\n"%(tstep))
f.write("variable            run_time equal %.1f                            # simulation run time [seconds]\n"%(runtime))
f.write("variable            frame_rate equal %.1f                          # dumping interval [seconds]\n"%(frate))
f.write("variable            seed equal %d                                  # random number generator seed\n"%(seed))
f.write('''
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
angle_coeff         * ${Kbend} 180.0

pair_style          hybrid/overlay zero 1.50 cosine/squared 1.50
pair_coeff          * * cosine/squared 0.00 1.00 1.10
pair_coeff          * * zero 1.50
pair_coeff          1 1 cosine/squared 1.00 1.00 1.00 wca
pair_coeff          1 2 cosine/squared 1.00 1.00 1.00 wca
pair_coeff          1 3 cosine/squared 1.00 1.00 1.00 wca
pair_coeff          2 2 cosine/squared 1.00 1.00 1.00 wca
pair_coeff          2 3 cosine/squared 1.00 1.00 1.00 wca
pair_coeff          3 3 cosine/squared 1.00 1.00 1.00 wca

# Nucleation reactions molecular templates
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
                react Nucleation all ${rstep} 0.900000 1.100000 mPreNucleation mPostNucleation Reactions/map_Nucleation.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create overlap 0.9 modify_create nuc yes         &
                react DimerOn all ${rstep} 0.900000 1.100000 mPreDimerOn mPostDimerOn Reactions/map_DimerOn.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create fit 1 modify_create overlap 0.9         &
                react TrimerOn all ${rstep} 0.900000 1.100000 mPreTrimerOn mPostTrimerOn Reactions/map_TrimerOn.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create fit 1 modify_create overlap 0.9         &
                react OligomerOn all ${rstep} 0.900000 1.100000 mPreOligomerOn mPostOligomerOn Reactions/map_OligomerOn.txt prob 1.0 ${seed} stabilize_steps ${stab_steps} modify_create fit 1 modify_create overlap 0.9         &
                react OligomerOff all ${rstep} 0.900000 1.100000 mPreOligomerOff mPostOligomerOff Reactions/map_OligomerOff.txt prob 1.0 ${seed} stabilize_steps ${stab_steps}         &
                react QuartomerOff all ${rstep} 0.900000 1.100000 mPreQuartomerOff mPostQuartomerOff Reactions/map_QuartomerOff.txt prob 1.0 ${seed} stabilize_steps ${stab_steps}         &
                react TrimerOff all ${rstep} 0.900000 1.100000 mPreTrimerOff mPostTrimerOff Reactions/map_TrimerOff.txt prob 1.0 ${seed} stabilize_steps ${stab_steps}         &
                react DimerOff all ${rstep} 0.900000 1.100000 mPreDimerOff mPostDimerOff Reactions/map_DimerOff.txt prob 1.0 ${seed} stabilize_steps ${stab_steps}

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

dump                1 all custom ${dump_time} output.xyz id mol type x y z v_vTailsTime v_vTailsP
dump_modify         1 format line "%d %d %d %.2f %.2f %.2f %.1f %.2f"

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
