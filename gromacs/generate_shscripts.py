import glob
import subprocess
import numpy as np

folder1 = glob.glob('./system_gRNA*')
folder2 = glob.glob('./system_cgI*')
folders = folder1 + folder2

structures = ['gRNA_loading.pdb','gRNA_ns.pdb','gRNA_nexus_del1.pdb',
              'gRNA_handle_del2.pdb','gRNA_handle_del1.pdb',
              'gRNA_nexus_del2.pdb','gRNA_nexus_mut1.pdb',
              'gRNA_nexus_xt1.pdb','gRNA_ws.pdb','cgIns_loading.pdb',
              'cgIws_loading.pdb','cgI_denovo_loading.pdb']

print(np.array(folders))
print(np.array(structures))

# create the systems and link the relevant structures
for i in range(0,len(structures)):
    structure = structures[i]
    folder = folders[i]
    cmds = ['ln -sf ../structures/'+structure+' '+folder+'/gRNA.pdb',
           'ln -sf ../structures/dCas9_loading.pdb '+folder+'/dCas9.pdb']
    for cmd in cmds:
        retval = subprocess.call(cmd,shell = True)

# Write bash script for doing equilibration
f = open('setupMD.sh','w')
f.write('#/bin/sh \n')
for i in range(0,len(structures)):
    folder = folders[i]
    cmd = 'python3 setup_system.py '+folder[1:]+'\n'
    f.write(cmd)
f.close()

# Write bash script for mdrun
f = open('runMD.sh','w')
f.write('#/bin/sh \n')
for i in range(0,len(structures)):
    folder = folders[i]
    cmd = 'gmx mdrun -v -deffnm '+folder+'/prod/prod \n'
    f.write(cmd)
f.close()

def backupData():
    files = glob.glob('./system*/prod/prod_nowater.xtc')
    for f in files:
        fname = f.split('/')[1][7:]
        cmd = 'mv '+f+' ./data/'+fname+f[-4:]
        print(cmd)
        subprocess.call(cmd, shell = True)

    files = glob.glob('./system*/prod/prod_nowater.pdb')
    for f in files:
        fname = f.split('/')[1][7:]
        cmd = 'mv '+f+' ./data/'+fname+f[-4:]
        print(cmd)
        subprocess.call(cmd, shell = True)


