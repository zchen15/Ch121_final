# zchen@caltech.edu
# 2019-04-08
# Some useful functions from rosetta and gromacs

import subprocess
import os
import glob
import time

import pandas # used for exporting PDB info
import Bio.PDB # for reading PDB files
import mdtraj # used for trajectory analysis
import numpy

import nglview

#############################################################################################################
# Useful shell and input output functions

# Wrapper function to submit shell commands and get its output
def submitShell(cmd):
    retval = subprocess.run(cmd, shell = True, stdout = subprocess.PIPE)
    print(retval.args)
    print(retval.stdout.decode('utf-8'))
    print(retval.returncode)

# Submit shell commands to background, such as doing a production run
def submitBackground(cmd):
    pid = subprocess.Popen([cmd], shell = True)
    print(cmd)
    print('Submitting job as pid ', pid)

# Rename residues in PDB to make it readable by the force field file
def readPDB(infile):
    ifile = open(infile, 'r')
    text = ifile.read()
    ifile.close()
    text = text.split('\n')
    data = []
    for line in text:
        x = [None] * 16
        if line[0:4]=='ATOM': # Record only ATOM information
            x[0] = 'ATOM'
            x[1] = int(line[6:11])
            x[2] = line[12:16].replace(' ','')
            x[3] = line[16]
            x[4] = line[17:20].replace(' ','')
            x[5] = line[21]
            x[6] = int(line[22:26])
            x[7] = line[26]
            x[8] = float(line[30:38])
            x[9] = float(line[38:46])
            x[10] = float(line[46:54])
            x[11] = float(line[54:60])
            x[12] = float(line[60:66])
            x[13] = line[72:76]
            x[14] = line[76:78].replace(' ','')
            x[15] = line[78:80]
            data.append(x)
    return pandas.DataFrame(data, columns=['ATOM','ATOM number','ATOM name','Alt loc',
                                       'Residue','Chain','Residue number','Residue insert code',
                                       'X','Y','Z','occupancy','temperature factor',
                                       'segment identifier','element','charge'])

# Write a table in pdb format
def writePDB(data, ofile):
    ofile = open(ofile,'w')
    for x in data:
        # Fucking bullshit PDB format which is not tab delimited
        fmt = '{0:6}{1:>5} {2:4}{3:1}{4:>3} {5:1}{6:>4}{7:1}   {8:>8}{9:>8}{10:>8}{11:>6}{12:>6}      {13:2}{14:>2}{15:>2}\n'.format(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15])
        ofile.write(fmt)
    ofile.close()

#############################################################################################################
# Useful rosetta functions

# Setup environmental variables so rosetta python scripts work
def setup_RosettaEnv(path='/opt/rosetta'):
    os.environ['ROSETTA'] = path
    os.environ['RNA_TOOLS'] = os.environ['ROSETTA']+'/tools/rna_tools'
    os.environ['PATH'] = os.environ['PATH']+':'+os.environ['RNA_TOOLS']+'/bin'
    os.environ['PYTHONPATH'] = os.environ['RNA_TOOLS']+'/bin'

def call_RNA_denovo(seq, dot, N, outfile, config=''):
    '''  
    This function calls the rna_denovo (FARFAR) procedure in rosetta.
    Avoid running this on denovo sequences longer than 100nt.
    Results are usually crap.
    
    seq = target sequence to build
    dot = secondary structure of target sequence
    N = number of structures to build
    outfile = name out silent output file
    configs = additional configurations to submit to rna_denovo binary
    
    Example:
    seq = 'aaacccttt'
    dot = '(((...)))'
    N = 10
    outfile = 'test.pdb'
    call_RNA_denovo(seq, dot, N, outfile)
    
    
    Other relevant configs
    -in:native                      Native PDB filename. Used for RMSD scoring calculations. [File]

    -minimize_rna                   High resolution optimize RNA after fragment assembly.[Boolean]

    -vary_geometry                  Vary bond lengths and angles (with harmonic constraints near 
                                    Rosetta ideal) for backbone and sugar degrees of freedom [Boolean]

    -cycles                         Number of Monte Carlo cycles.[default 10000]. [Integer]

    -bps_moves                      Base pair step moves. For adjacent base pairs within stems or that 
                                    are obligate pairs, draw sequence-matched fragments that encompass
                                    both pairs. Adjacent means that base pairs have contiguous residues
                                    on one strand, and at most 3 intervening residues on the other.

    -output_lores_silent_file       If high resolution minimizing, output intermediate low resolution models.
                                    [Boolean]

    -dump                           Generate pdb output. [Boolean]
    
    -vall_torsions                  Source of RNA fragments. [default: 1jj2.torsions]. [Boolean]
    
    -jump_library_file              Source of base-pair rigid body transformations if base pairs are specified.
                                    [default: 1jj2_RNA_jump_library.dat] [String]

    -obligate_pair                  Residue pairs that must form a base pair (possibly non canonical)
    
    -secstruct_general              Specification of -obligate_pair in dot-parens format

    -obligate_pair_explicit         Residue pairs that must form a base pair, with  specification of base 
                                    edges (W/H/S/X) and orientation (A/P/X for antiparallel/
                                    parallel/unknown; C/T/X allowed too for cis/trans)

    -cst_file                       Specify constraints (typically atom pairs) in Rosetta-style constraint
                                    file. [String]

    -output_lores_silent_file       if doing full-atom minimize, also save models after fragment assembly
                                    but before refinement (file will be called *.LORES.out) [Boolean]
    
    -dump                           output pdbs that occur during the run, even if using silent file output.
    
    Advanced Configs
    -s                              Input PDBs to be used as fixed 'chunks' in fragment assembly

    -in:file:silent                 List of input files (in 'silent' format) that specify potential template
                                    structures or 'chunks'

    -input_res                      Positions at which 'chunks' are applied. If there is more than one chunk
                                    file, specify indices for the first file and then the second file, etc.
                                    (Used to be called -chunk_res.)

    -fixed_stems                    Seed each stem with a Watson-Crick base pair instead of having the 
                                    strands find each other

    -cutpoint_closed                Positions at which to force transient chainbreaks (may be needed
                                    if you get fold-tree errors)

    -cutpoint_open                  Positions at which strands end (better to specify separate strands
                                    in FASTA file, or with spaces between strings in sequence)

    -data_file                      RDAT or legacy-format file with RNA chemical mapping data
    
    Less Common Configs
    -filter_lores_base_pairs        Filter for models that satisfy structure parameters.
                                    [Boolean] True by default.

    -in:database                    Path to rosetta databases. Default is based on location of 
                                    rosetta executables. [PathVector]
    
    -output_res_num                 Numbering (and chain) of residues in output PDB or silent file. 
                                    Better to specify in headers in .fasta file.
    
    -staged_constraints             Apply constraints in stages depending on sequence separation
    
    -close_loops                    Attempt closure across chainbreaks by cyclic coordinate descent
                                    after fragment moves [Boolean] Defaults to true.
    '''
    start = time.time()
    cmd = ['./bin/rna_denovo.linuxgccrelease '+config+' -sequence "'+ seq.lower() +'" -secstruct "'+ dot 
           +'" -nstruct '+ str(N) +' -out:file:silent '+ outfile+' -tag '+outfile]
    print(cmd)
    retval = subprocess.call(cmd, shell = True)
    print('retval = ',retval,' elapse = ',time.time()-start)

def call_RNA_stepwiseMC(pdb, target, fastafile, outfile, configs = ''):
    '''
    This function calls the stepwise monte carlo procedure in rosetta, which is used to build
    protein and RNA loops of unknown structure. If you know the approximate base pairing structure,
    use call_RNA_denovo() instead.
    
    Inputs:
    pdb = seed helices
    target = sequence to build of which some must match the pdb file
    fastafile = file to write target sequence info in fasta format
    outfile = name of rosetta silent file as output
    
    Returns:
    The function generates a rosetta silent file which contains the optimized structures.
    Use call_RNA_getPDB() to extract the pdb structure
    
    Example:
    
    call_RNA_stepwiseMC('test.pdb', 'aaaacccc', 'seq.fa', 'stepwise.out')
    
    Other relevant configs
    -stepwise:monte_carlo:cycles                    Species how many monte carlo cycles to do
                                                    [Default = 50] 

    -stepwise:protein:allow_virtual_side_chains     Allow virtual side chains in packing
    
    -stepwise:rna:erraser                           Use KIC sampling for loop modeling
    
    -align_pdb                                      pdb file to align to
    
    -add_delete_frequency                           Frequency of add/delete vs resampling
                                                    [Default = 0.5]

    -minimize_single_res_frequency                  Frequency to minimize current residue before global minimization
                                                    [Default = 0]
    
    -temperature                                    Sets the monte carlo temperature
                                                    [Default = 1]
    
    '''
    start = time.time()
    # Write alignment info to fasta file
    print('Writing sequence to fasta file')
    f = open(fastafile,'w')
    f.write('> TargetSequence\n')
    f.write(target.lower()+'\n')
    f.close()
    
    # execute the stepwise function
    cmd = ['./bin/stepwise.linuxgccrelease '+configs+' -fasta '+fastafile+' -s '+pdb+' -out:file:silent '+outfile
           +' -score:rna_res_level_energy7beta']
    print(cmd)
    retval = subprocess.call(cmd, shell = True)
    print('retval = ',retval,' elapse = ',time.time()-start)

    if retval!=0: # catch errors
        print('Stepwise function failed')

def call_RNA_helix(seq, outfile, configs='-minimize_rna'):
    '''    
    This function builds an RNA helix based on given input sequences.
    This function takes two RNA sequences separated by a space. These sequences should be the same length.
    
    Example:
    call_RNA_helix('aaaa tttt', 'test.pdb')

    '''
    start = time.time()
    cmd = './bin/rna_helix.linuxgccrelease '+configs+' -seq '+seq.lower()+' -o '+outfile
    print(cmd)
    retval = subprocess.call(cmd, shell = True)
    print('retval = ',retval,' elapse = ',time.time()-start)

def call_RNA_thread(target, template, pdb, fastafile, outfile, resN = 0):
    '''    
    This function fits a target sequence to a given template pdb.
    
    This function can also be used to mutate specific residues on
    the template pdb. However, the target and template sequences must be same lengths
    
    target = aligned target sequence
    template = aligned template sequence
    NOTE: target and template sequences must be aligned to each other
    
    pdb = template pdb file
    fastafile = file to write fasta alignment info
    outfile = output pdb file
    resN = offset numbering for output pdb
    
    Example:
    
    pdb = 'gRNA_loading.pdb'
    seq1 = 'aacuuucaguuuagcggucuguuuuagagcuagaaauagcaaguuaaaauaaggcuaguccgACACAAGGGG----AAAUUAACAACACAACACACACAACACAGGccccgg'
    seq2 = '----------GAUGAGACGCguuuuagagcuagaaauagcaaguuaaaauaaggcuaguccguuaucaacuugaaa------------------------------aagugu'
    call_RNA_thread(seq1,seq2,pdb,'aligned.fa','cgI.pdb')

    '''
    start = time.time()
    # Write alignment info to fasta file
    f = open(fastafile,'w')
    f.write('> TargetSequence\n')
    f.write(target.lower()+'\n')
    f.write('> TemplateSequence\n')
    f.write(template.lower()+'\n')
    f.close()
    
    # calls the rosetta binary
    cmd = ['./bin/rna_thread.linuxgccrelease -in:file:fasta '+fastafile
           +' -s '+pdb+' -o '+outfile+' -seq_offset '+str(resN)]
    print(cmd)
    retval = subprocess.call(cmd, shell = True)
    print('retval = ',retval,' elapse = ',time.time()-start)
    
    if retval==0: # Catch errors
        return 0
    else: # Remove fasta file
        print('rna_thread failed')
    
# put together RNA fragments
def call_RNA_graft(frags, outfile):
    '''
    This function grafts RNA fragments into a full RNA structure.
    
    Example:
    
    call_RNA_graft(['frag1.pdb','frag2.pdb'],'test.pdb')
    
    '''
    start = time.time()
    if isinstance(frags,str): # Catch stupid errors like giving only 1 fragment
        print('only one fragment present')
        cmd = './bin/rna_graft.linuxgccrelease -s '+frags+' -o '+outfile
    else:
        cmd = './bin/rna_graft.linuxgccrelease -s '
        for f in frags:
            cmd = cmd + ' ' + f
        cmd = cmd + ' -o '+ outfile
    print(cmd)
    retval = subprocess.call(cmd, shell = True)
    print('retval = ',retval,' elapse = ',time.time()-start)

def call_RNA_splice(infile, resnum, outfile):
    '''
    This function calls a rosetta script which splices specific residues out of a pdb file.
    
    Example:
    call_RNA_splice('input.pdb','A:1-10 B:20-25','output.pdb')
    
    This takes input.pdb and splices out chain A residues 1-10 and chain B residues 20-25.
    The spliced out structures are saved in output.pdb
    
    '''
    cmd = './rna_tools/bin/pdbslice.py '+infile+' -subset '+resnum+' subset_'
    print(cmd)
    retval = subprocess.run(cmd, shell = True)
    
    # Check for errors on execution
    if retval==0:
        subprocess.run('mv subset_'+infile+' '+outfile, shell = True)
    else: # Moves the output to outfile
        print('error on pdbslice.py')

# Splice out a chain from a pdb file
def call_RNA_getchain(infile, chain, outfile):
    '''
    This function calls a rosetta script which extracts a specific chain from a pdb file.
    
    Example:
    call_RNA_getchain('test.pdb','A','out.pdb')
    
    '''
    cmd = 'python2 ./rna_tools/bin/extract_chain.py '+infile+' '+chain
    print(cmd)
    subprocess.run(cmd, shell = True)
    
    # Rename the script output to the user desired output file
    text = infile.split('.')
    out = text[0]+chain+'.'+text[1]
    subprocess.run('mv '+out+' '+outfile, shell = True)
    
def call_RenumberPDB(infile, resN):
    '''
    This function calls renumber_pdb_in_place.py to renumber the residues of a pdb file in place.
    This function is typically used in combination with RNA grafting because rna_graft puts together
    fragments based on the residue locations.
    
    Example:
    call_RenumberPDB('test.pdb','A:1-20 B:1-20')
    
    '''
    cmd = './rna_tools/bin/renumber_pdb_in_place.py '+infile+' '+resN
    print(cmd)
    result = subprocess.run(cmd, shell = True)
    print('retval = ',result)

def call_getSequence(file):
    '''    
    This function calls a rosetta script which parses a pdb file and
    return all residues in the pdb file.
    '''
    cmd = './rna_tools/bin/get_sequence.py '+file
    print(cmd)
    result = subprocess.run(cmd, shell = True, stdout = subprocess.PIPE)
    if result.returncode==0:
        return str(result.stdout)[2:-3]
    else:
        print('Getting sequence execution failed')
        return -1

def call_RNA_getPDB(infile, outfile):
    '''    
    This function extracts pdb files from rosetta silent file format.
    The function returns the list of these files with the output prefix.
    
    Example:
    call_getPDB('test.out', 'test')
    
    Returns = ['testS_00001.pdb','testS_00002.pdb','testS_00003.pdb',...]
    '''
    start = time.time()
    cmd = ['./bin/extract_pdbs.linuxgccrelease -in:file:silent '+infile
           +' -out:prefix '+outfile+' -in:file:silent_struct_type rna']
    print(cmd)
    retval = subprocess.call(cmd, shell = True)
    print('retval = ',retval,' elapse = ',time.time()-start)
    
    if retval==0: # catch errors
        return glob.glob(outfile+'*.pdb')
    else: # if it works return the pdb file list
        print('Extracting pdb failed')
        return -1

#############################################################################################################
# Useful gromacs functions

def gromacs_recenter_structure(infile, topo, outfile, select = 'Protein System', configs = '-center -pbc mol -ur compact'):
    '''
    This function calls gromacs trjconv to recenter the pdb or trajectory files. If structures
    drift during simulation, they can hit the periodic boundary position and look really weird.
    Recentering the structure will make it look better.

    Inputs:
    infile = input structure or topology
    topo = topology of the input structure file
    outifle = name of output file

    Example:
    recenter_structure('input.pdb','input.top','output.pdb')

    '''
    cmd = 'echo '+select+' | ' # recenter on the protein group
    cmd = cmd + 'gmx trjconv '+configs+' -s '+topo+' -f '+infile+' -o '+outfile
    submitShell(cmd)
    
def gromacs_convert_to_PDB(infile, topo, outfile, select = 'non-Water', configs = ''):
    '''
    This function calls gromacs trjconv to recenter the pdb or trajectory files. If structures drift
    during simulation, they can hit the periodic boundary position and look really weird. Recentering
    the structure will make it look better.

    Inputs:
    infile = input structure or topology
    topo = topology of the input structure file
    outifle = name of output file

    Example:
    recenter_structure('input.pdb','input.top','output.pdb')

    '''
    cmd = 'echo '+select+' | ' # select non-waters and write to new file
    cmd = cmd + 'gmx trjconv '+configs+' -s '+topo+' -f '+infile+' -o '+outfile
    submitShell(cmd)

#############################################################################################################
# Useful nglview functions
def add_structure(file, view=None):
    '''
    Takes pdb file as input and adds it to a new nglview widget
    if a view widget is not given, then instantiate one
    '''
    S1 = Bio.PDB.PDBParser().get_structure(file, file)
    nS1 = nglview.BiopythonStructure(S1)
    # if view is empty, then make one
    if view == None:
        view = nglview.widget.NGLWidget(gui=True)
    view.add_structure(nS1)
    return view

def color_structure(view, rep_type, sele, resN, color):
    '''
    This function colors a structure in an nglview object

    Inputs:
    view = nglview object
    rep_type = cartoon, ball+stick, or etc...
    sele = protein, nucleic, not protein, and etc
    resN = vector of integers representing boundaries between domains
    color = colors for each domain

    Example:
    gRNA_domainN = [1,10,23,26,40,53,59,63,93,98]
    gRNA_domainColor = ['red','blue','white','blue','orange','purple','green','white','green']
    color_structure(view, 'cartoon', 'nucleic', gRNA_domainN, gRNA_domainColor)
    '''
    for i in range(0,len(color)):
        view.add_representation(rep_type, sele+' and '+str(resN[i])+'-'+str(resN[i+1]), color = color[i])

# generate list of hydrogen bonds from mdtraj information
def mdtraj_get_HBondList(traj):
    '''
    
    '''
    # compute hbonds using wernet nilsson
    hbond = mdtraj.wernet_nilsson(traj, periodic=False)
    data = []
    # get dataframe from traj and use it to compute chain id
    table, bonds = traj.top.to_dataframe()
    
    for frame in range(0,len(hbond)):
        print('Computing frame ', frame)
        for bond in range(0,len(hbond[frame])):
            # get atom index
            d_i = hbond[frame][bond][1]
            a_i = hbond[frame][bond][2]
            
            # get donor and acceptor names
            donor = traj.topology.atom(d_i)
            acceptor = traj.topology.atom(a_i)

            # get positions
            d_pos = traj.xyz[frame][d_i]
            a_pos = traj.xyz[frame][a_i]
            
            # get chain location
            d_chain = int(table[table['serial']==d_i]['chainID'])
            a_chain = int(table[table['serial']==a_i]['chainID'])
            
            # get donor residue numbers
            d_resn = int(table[table['serial']==d_i]['resSeq'])
            a_resn = int(table[table['serial']==a_i]['resSeq'])
                      
            # save the data
            data.append([frame, d_chain, d_resn, donor, d_pos, a_chain, a_resn, acceptor, a_pos])
    df = pandas.DataFrame(data, columns = ['frame','donor chain','donor resnum','donor','donor position (nm)',
                                       'acceptor chain','acceptor resnum','acceptor','acceptor position (nm)'])
    # compute the bond length
    x = df['donor position (nm)'] - df['acceptor position (nm)']
    x = numpy.array([i for i in x.values]) # convert numpy array                                  
    df['length (nm)'] = numpy.linalg.norm(x, axis=1)
    return df        

def drawHbond(view, data):
    '''
    This function takes an nglview object and hydrogen bond pandas dataframe
    as input and draws arrows representating hydrogen bonds for each frame in
    the nglview object.
    
    hydrogen bond dataframe is generated by mdtraj_get_HBondList
    '''
    frames = numpy.unique(data['frame'])
    for f in frames:
        view.frame = int(f)
        df = data[data['frame']==f]
        for i in range(0,len(df)):
            # draws arrows for hydrogen bonds
            x = df.iloc[i]['donor position (nm)']*10 # scale by 10 because nglview draws in angstroms
            pos1 = [x[0], x[1], x[2]]
            x = df.iloc[i]['acceptor position (nm)']*10 # scale by 10 because nglview draws in angstroms
            pos2 = [x[0], x[1], x[2]]
            color = [0,255,0] # make it green
            label = str(df.iloc[i]['donor']) + '-' + str(df.iloc[i]['acceptor'])
            view.shape.add_arrow(pos1, pos2, color , 0.1, label)
            # draw the sidechains as licorice
            resn1 = df.iloc[i]['donor resnum']
            resn2 = df.iloc[i]['acceptor resnum']
            
            view.add_representation('licorice',selection='protein and '+str(resn1))
            view.add_representation('licorice',selection='protein and '+str(resn2))
        
def add_axis(view):
    view.shape.add_arrow([0,0,0],[0,0,100],[0,0,255],5.0,'Z')
    view.shape.add_arrow([0,0,0],[0,100,0],[0,255,0],5.0,'Y')
    view.shape.add_arrow([0,0,0],[100,0,0],[255,0,0],5.0,'X')