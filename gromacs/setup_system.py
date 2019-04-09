
# coding: utf-8

# In[ ]:


# General tools
import numpy as np
import scipy as sp
import pandas as pd

import parmed.gromacs

import os
import sys

from mdtools import *


# In[2]:

cwd = os.getcwd()+'/'
newcwd = cwd+sys.argv[1]
os.chdir(newcwd)
print(os.getcwd())

# In[ ]:


# Generating topology file for gRNA and tip3p water
RNA = {'U':'RU', 'G':'RG', 'C':'RC', 'A':'RA'} # need to rename RNA residues from rosetta style so gromacs forcefields can understand it
df = readPDB('gRNA.pdb')
x = df['Residue'].values # Rename residues for amber compatibility
for i in range(0,len(x)):
    if x[i] in RNA:
        x[i] = RNA[x[i]]
df['Chain'] = 'B' # rename to chain B
df = df[df['Residue number']>1] # Toss first residue
df = df[3:] # Cut off at O5' due to some rosetta issues
writePDB(df.values, 'gRNA_fixed.pdb')

prefix = 'gRNA'
# Generate the topology file
cmd = ['gmx pdb2gmx -f '+prefix+'_fixed.pdb -ff amber99sb-ildn -ignh -o '
       +prefix+'.gro -p '+prefix+'.top -i '+prefix+'_posre.itp -water tip3p']
submitShell(cmd)

prefix = 'dCas9'
# Generate the topology file
cmd = ['gmx pdb2gmx -f '+prefix+'.pdb -ff amber99sb-ildn -ignh -o '
       +prefix+'.gro -p '+prefix+'.top -i '+prefix+'_posre.itp -water tip3p']
submitShell(cmd)


# In[ ]:


# Merge the two topology files
def readTopFile(file):
    f = open(file,'r')
    out = f.read().split('\n')
    f.close()
    return out

def writeTopFile(lines, file):
    f = open(file,'w')
    for line in lines:
        f.write(line+'\n')
    f.close()

# Read the topology files
gRNA_top = readTopFile('gRNA.top')
dCas9_top = readTopFile('dCas9.top')
# Truncate and write itp files
if os.path.exists('./topology') == False:
    subprocess.call('mkdir topology', shell = True)
else:
    subprocess.call('rm ./topology/*', shell = True)

subprocess.call('mv gRNA_posre.itp ./topology/', shell = True)
subprocess.call('mv dCas9_posre.itp ./topology/', shell = True)

writeTopFile(gRNA_top[21:-27], './topology/gRNA.itp')
writeTopFile(dCas9_top[21:-27], './topology/dCas9.itp')
# add molecules together
molecules = ''
molecules = molecules + dCas9_top[-2] + '\n'
molecules = molecules + gRNA_top[-2] + '\n'

# include optional position restraint for gRNA, dCas9, and water
system = '''
; Include forcefield parameters
#include "amber99sb-ildn.ff/forcefield.itp"

; Include dCas9 topology
#include "./topology/dCas9.itp"'

; Include Position restraint file
#ifdef POSRES_DCAS9
#include "./topology/dCas9_posre.itp"
#endif

; Include gRNA topology
#include "./topology/gRNA.itp"'

#ifdef POSRES_GRNA
#include "./topology/gRNA_posre.itp"
#endif

; Include generic topology for ions
#include "amber99sb-ildn.ff/ions.itp"

; Include water topology
#include "amber99sb-ildn.ff/tip3p.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

[ system ]
; Name
dCas9_gRNA

[ molecules ]
; Compound        #mols
'''

# Write the new system topology file
f = open('system.top','w')
f.write(system)
f.write(molecules)
f.close()

# Combine the two gromacs files
dCas9_gro = parmed.gromacs.GromacsGroFile.parse('dCas9.gro')
gRNA_gro = parmed.gromacs.GromacsGroFile.parse('gRNA.gro')
parmed.gromacs.GromacsGroFile.write(dCas9_gro + gRNA_gro,'system.gro')

# Remove the old topology and gro files
clear_files = True
if clear_files:
    subprocess.call('rm gRNA.top', shell = True)
    subprocess.call('rm dCas9.top', shell = True)
    subprocess.call('rm gRNA.gro', shell = True)
    subprocess.call('rm dCas9.gro', shell = True)


# In[ ]:


# Define the simulation box
cmd = 'gmx editconf -f system.gro -o system.gro -bt dodecahedron -d 1.0'
#cmd = 'gmx editconf -f system.gro -o system.gro -c -d 1.0 -bt cubic'
submitShell(cmd)

# solvate and add ions
# 10mM NaCl for bacteria
# 100mM for blood plasma
# http://book.bionumbers.org/what-are-the-concentrations-of-different-ions-in-cells/
submitShell('cp system.top solvated.top')
cmd = 'gmx solvate -cp system.gro -cs spc216 -p solvated.top -o solvated.gro'
submitShell(cmd)

protocol = '''
; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
coulombtype	    = cutoff	; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; long range electrostatic cut-off
rvdw		    = 1.0		; long range Van der Waals cut-off
pbc             = xyz 		; Periodic Boundary Conditions
'''
f = open('ions.mdp','w')
f.write(protocol)
f.close()

cmd = 'gmx grompp -f ions.mdp -c solvated.gro -p solvated.top -o ions.tpr'
submitShell(cmd)

cmd = 'echo SOL |'
cmd = cmd + 'gmx genion -s ions.tpr -o solvated.gro -p solvated.top -pname NA -nname CL -neutral -conc 0.01'
submitShell(cmd)

# Remove the old topology and gro files
if clear_files:
    subprocess.call('rm ions.tpr', shell = True)
    subprocess.call('rm ions.mdp', shell = True)


# In[ ]:


# Setup minimization
protocol = '''
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		        ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		    ; Method to determine neighbor list (simple, grid)
rlist		    = 1.2		    ; Cut-off for making neighbor list (short range forces)
coulombtype	    = PME		    ; Treatment of long range electrostatic interactions
rcoulomb	    = 1.2		    ; long range electrostatic cut-off
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw		    = 1.2		    ; long range Van der Waals cut-off
pbc             = xyz 		    ; Periodic Boundary Conditions
DispCorr        = no
'''

# Clear the minimization folder
if os.path.exists('./min') == False:
    subprocess.call('mkdir min', shell = True)
else:
    subprocess.call('rm ./min/*', shell = True)

f = open('./min/em.mdp','w')
f.write(protocol)
f.close()

cmd = 'gmx grompp -f ./min/em.mdp -c solvated.gro -p solvated.top -o ./min/em.tpr -maxwarn 10'
submitShell(cmd)
print('Running energy minimization')
submitShell('gmx mdrun -v -deffnm ./min/em')


# In[3]:


# Setup NVT
protocol = '''
title                   = Equilibration NVT 
define                  = -DPOSRES_DCAS9 -DPOSRES_GRNA ; position restrain the protein and RNA

; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 5000     ; 2 * 5000 = 10 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 0         ; save coordinates every 1.0 ps
nstvout                 = 0         ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps

; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme

; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 310     310           ; reference temperature, one for each group, in K

; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 310       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed

'''
# Clear the folder
if os.path.exists('./nvt') == False:
    subprocess.call('mkdir nvt', shell = True)
else:
    subprocess.call('rm ./nvt/*', shell = True)

f = open('./nvt/nvt.mdp','w')
f.write(protocol)
f.close()

cmd = 'gmx grompp -f ./nvt/nvt.mdp -c ./min/em.gro -r ./min/em.gro -p solvated.top -o ./nvt/nvt.tpr -maxwarn 10'
submitShell(cmd)
print('Running nvt')
submitShell('gmx mdrun -v -deffnm ./nvt/nvt')


# In[4]:


# Setup NPT
protocol = '''
title                   = equilibration NPT 
;define                  = -DPOSRES_GRNA -DPOSRES_DCAS9 ; position restrain the protein

; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 5000     ; 2 * 5000 = 10 ps
dt                      = 0.002     ; 2 fs

; Output control
nstxout                 = 0         ; save coordinates every 1.0 ps
nstvout                 = 0         ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps

; Bond parameters
continuation            = yes       ; Restarting after NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme

; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 310     310           ; reference temperature, one for each group, in K

; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com

; Periodic boundary conditions
pbc                     = xyz  ; 3-D PBC

; Velocity generation
gen_vel                 = no        ; Velocity generation is off 

'''
# Clear the folder
if os.path.exists('./npt') == False:
    subprocess.call('mkdir npt', shell = True)
else:
    subprocess.call('rm ./npt/*', shell = True)

f = open('./npt/npt.mdp','w')
f.write(protocol)
f.close()

cmd = 'gmx grompp -f ./npt/npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro -p solvated.top -o ./npt/npt.tpr -maxwarn 10'
submitShell(cmd)
print('Running npt')
submitShell('gmx mdrun -v -deffnm ./npt/npt')


# In[5]:


# Setup production MD run
protocol = '''
title                   = Production MD run
; define                  =  DPOSRES_GRNA -DPOSRES_DCAS9 ; position restrain the protein

; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 5000000    ; 2 * 5000000 = 1000 ps (10 ns)
dt                      = 0.002     ; 2 fs

; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying 
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps
compressed-x-grps       = System    ; save the whole system

; Bond parameters
continuation            = yes       ; Restarting after NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

; Neighborsearching
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)

; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 310     310           ; reference temperature, one for each group, in K

; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1

; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC

; Dispersion correction
DispCorr                = EnerPres  ; account for cut-off vdW scheme

; Velocity generation
gen_vel                 = no        ; Velocity generation is off 

'''
# Clear the folder
if os.path.exists('./prod') == False:
    subprocess.call('mkdir prod', shell = True)
else:
    subprocess.call('rm ./prod/*', shell = True)

f = open('./prod/prod.mdp','w')
f.write(protocol)
f.close()

cmd = 'gmx grompp -f ./prod/prod.mdp -c ./npt/npt.gro -p solvated.top -o ./prod/prod.tpr -maxwarn 10'
submitShell(cmd)
#print('Running production run')
#submitBackground('gmx mdrun -v -deffnm ./prod/prod')


# In[ ]:


# Setup simulated annealing run
protocol = '''
title		= Simulated Annealing 

;define		= -DPOSRES_GRNA -DPOSRES_DCAS9	; restraints

; Run parameters
integrator	= md		; leap-frog integrator
dt			= 0.002		; 2 fs
nsteps		= 50000	; 100 ps

; Output control
nstxout					= 1000		; save coordinates every 2 ps
nstvout					= 1000 		; save velocities every 2 ps
nstfout					= 1000		; save forces every 2 ps
nstenergy				= 5000		; save energies every 2 ps
nstlog					= 5000      ; update log file every 10.0 ps
nstxout-compressed		= 5000      ; save compressed coordinates every 10.0 ps
compressed-x-grps		= System    ; save the whole system

; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = all-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

; Neighborsearching
nstlist		= 5		    ; 10 fs
ns_type		= grid 		; search neighboring grid cells
rlist		= 1.2		; short-range neighborlist cutoff (nm)
rcoulomb	= 1.2		; short-range electrostatic cutoff (nm)
rvdw		= 1.2		; short-range van der Waals cutoff (nm)

; Electrostatics
coulombtype		= PME		; Particle Mesh Ewald for long-range electrostatics
pme_order		= 4			; cubic interpolation
fourierspacing  = 0.16		; grid spacing for FFT

; Dispersion correction
DispCorr = EnerPres ; account for cut-off vdW scheme

; Temperature coupling is on in three groups
Tcoupl	 	= Berendsen					; Weak coupling
tc_grps		= System					; three coupling groups - more accurate
tau_t		= 0.1 						; time constant, in ps
ref_t		= 320 						; reference temperature, one for each group, in K

; Pressure coupling
Pcoupl				= Berendsen				; Weak coupling
Pcoupltype			= isotropic				; uniform scaling of x-y vectors, independent z
tau_p				= 0.5					; time constant, in ps
ref_p				= 1.0					; reference pressure, x-y, z (in bar)
compressibility		= 4.5e-5				; isothermal compressibility of water, bar^-1
refcoord_scaling	= com

; Periodic boundary conditions are on in all directions
pbc = xyz           ; 3-D PBC

; Simulated annealing
annealing         = single            ; single sequence of points for each T-coupling group
annealing_npoints = 2                 ; two points - start and end temperatures
annealing_time    = 0   500           ; time frame of heating - heat over period of 500 ps
annealing_temp    = 320 310           ; start and end temperatures

; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 320       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed

'''

'''
# Clear the minimization folder
if os.path.exists('./anneal') == False:
    subprocess.call('mkdir anneal', shell = True)
else:
    subprocess.call('rm ./anneal/*', shell = True)

f = open('./anneal/anneal.mdp','w')
f.write(protocol)
f.close()

cmd = 'gmx grompp -f ./anneal/anneal.mdp -c ./min/em.gro -p system.top -o ./anneal/anneal.tpr -maxwarn 10'
#submitShell(cmd)
print('Running annealing')
#submitShell(cmd)
'''


# In[ ]:
