#!/usr/bin/env python
import os
import subprocess

assert 0

burg = '/burg/berkelbach/users/eav2136'
moto = '/moto/berkelbach/users/eav2136'

material_components = {
        'c':['C'], 
        'si':['Si'], 
        'sic':['Si', 'C'], 
        'bn':['B', 'N'], 
        'bp':['B', 'P'], 
        'aln':['Al', 'N'], 
        'alp':['Al', 'P'], 
        'mgo':['Mg', 'O'], 
        'mgs':['Mg', 'S'], 
        'lih':['Li', 'H'], 
        'lif':['Li', 'F'], 
        'licl':['Li', 'Cl']}

def write_basis(formula, material_components, basis_file):
    basis = '{'
    for element_index, element in enumerate(material_components[formula]):
        if element_index > 0:
            basis = basis + ', '
        basis = basis + '\'' + element + '\': parse_nwchem.load(\'' + basis_file + '\', \'' + element + '\')'
    basis = basis + '}'
    return basis

# Write import headers block
def write_header(openfile):
    openfile.write('''#!/usr/bin/env python
import numpy as np
from pyscf.gto.basis import parse_nwchem
from pyscf.pbc import gto, scf, cc
from pyscf.pbc.cc import eom_kccsd_rhf
from pyscf.pbc.tools import pyscf_ase, lattice
from pyscf.pbc.cc.eom_kccsd_rhf import _IMDS
from pyscf.pbc.scf import chkfile
import os
import h5py

prefix = os.getcwd()
imds = None
    ''')

# Write Cell function block
def write_cell(openfile, formula, scaled_center, basis_file):
    basis = write_basis(formula, material_components, basis_file)
    openfile.write('''
##############################
# Create a "Cell"
##############################

cell = gto.Cell()
# Candidate formula of solid: c, si, sic, bn, bp, aln, alp, mgo, mgs, lih, lif, licl
formula = '{}'
scaled_center = {}
ase_atom = lattice.get_ase_atom(formula)
cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atom)
cell.a = ase_atom.cell
cell.unit = 'B'
cell.basis = {}
cell.pseudo = 'gthhfrev'
cell.verbose = 7
cell.build()
    '''.format(formula, scaled_center, basis))

# Write SCF block
def write_scf(openfile, kdensity, base_name):
    openfile.write('''
##############################
#  K-point SCF 
##############################

kdensity = {}
kmesh = [kdensity, kdensity, kdensity] 
kpts = cell.make_kpts(kmesh, scaled_center=scaled_center)
mymf = scf.KRHF(cell, kpts=kpts, exxdiv='ewald')
mymf = mymf.density_fit()

# save ERIs if not exist
cderi_h5name = '{}' + '_cderi.h5'
cderi_h5name = os.path.join(prefix, cderi_h5name)
mymf.with_df._cderi_to_save = cderi_h5name

if os.path.isfile(cderi_h5name):
    mymf.with_df._cderi = cderi_h5name

chkfile_name = '{}' + '.chk'
chkfile_name = os.path.join(prefix, chkfile_name)
mymf.chkfile = chkfile_name

if os.path.isfile(chkfile_name):
    chk_cell, chk_scf = chkfile.load_scf(chkfile_name)
    mymf.mo_coeff = chk_scf['mo_coeff']
    mymf.mo_energy = chk_scf['mo_energy']
    mymf.mo_occ = chk_scf['mo_occ']
    mymf.e_tot = chk_scf['e_tot']
else:
    ekrhf = mymf.kernel()
    '''.format(kdensity, base_name, base_name))

# Write CCSD block
def write_ccsd(openfile, base_name):
    openfile.write('''
##############################
# K-point CCSD
##############################

cc_h5name = '{}' + '_cc.h5'
cc_h5name = os.path.join(prefix, cc_h5name)

frozen = None
mycc = cc.KRCCSD(mymf, frozen=frozen)

if kdensity > 3:
    mycc.max_cycle = 2

if os.path.isfile(cc_h5name):
    print('Reading t1, t2 from {{}}'.format(cc_h5name))
    t1 = None
    t2 = None
    with h5py.File(cc_h5name, 'r') as fin:
        t1 = fin['t1'][:]
        t2 = fin['t2'][:]
    mycc.t1 = t1
    mycc.t2 = t2

ekrcc, t1, t2 = mycc.kernel()

if os.path.isfile(cc_h5name):
    print('Saving t1, t2 to {{}}'.format(cc_h5name))
    with h5py.File(cc_h5name, 'a') as fout:
        fout['t1'][:] = t1
        fout['t2'][:] = t2
else:
    print('Saving t1, t2 to {{}}'.format(cc_h5name))
    with h5py.File(cc_h5name, 'w') as fout:
        fout.create_dataset('t1', data=t1)
        fout.create_dataset('t2', data=t2)
    '''.format(base_name))

# Read CCSD t1 and t2 amplitudes
def read_ccsd(openfile, base_name):
    openfile.write('''
##############################
# K-point CCSD
##############################

cc_h5name = '{}' + '_cc.h5'
cc_h5name = os.path.join(prefix, cc_h5name)

frozen = None
mycc = cc.KRCCSD(mymf, frozen=frozen)

if os.path.isfile(cc_h5name):
    print('Reading t1, t2 from {{}}'.format(cc_h5name))
    t1 = None
    t2 = None
    with h5py.File(cc_h5name, 'r') as fin:
        t1 = fin['t1'][:]
        t2 = fin['t2'][:]
    mycc.t1 = t1
    mycc.t2 = t2
    '''.format(base_name))

# Write EOM intermediates

def write_imds_ip(openfile):
    openfile.write('''
if imds == None:
    imds = _IMDS(mycc)
    imds.make_ip()
else:
    imds.make_ip()
            ''')

def write_imds_ea(openfile):
    openfile.write('''
if imds == None:
    imds = _IMDS(mycc)
    imds.make_ea()
else:
    imds.make_ea()
            ''')

def save_imds(openfile):
    openfile.write('''
print('Saving imds to {}'.format(cc_h5name))
with h5py.File(cc_h5name, 'a') as fout:
    print("Saving imds to a group called 'imds'")
    grp = fout.create_group('imds')
    for k, v in imds.__dict__.items():
        if k[0] == "F" or k[0] == "W" or k[0] == "L":
            print("imds key to save:", k, ", type:", type(v))
            grp.create_dataset(k, data=v)
    ''')

def read_imds(openfile):
    openfile.write('''
# rebuild imds with necessary information
imds_dict = chkfile.load(cc_h5name, 'imds')
imds = _IMDS(mycc)
for k, v in imds_dict.items():
    setattr(imds, k, v)
imds_dict = None
    ''')

def write_eomip(openfile, nroots, base_name):
    openfile.write('''
##############################
# EOM-IP-KRCCSD
##############################

# number of roots requested
# index(indices) of targetted k point(s) 
kptlist = [0]

myeom = eom_kccsd_rhf.EOMIP(mycc)
eip, vip = myeom.ipccsd(nroots={}, imds=imds, kptlist=kptlist)

eomip_h5name = '{}' + '_eomip.h5'
eomip_h5name = os.path.join(prefix, eomip_h5name)

print('Saving eip, vip to {{}}'.format(eomip_h5name))
with h5py.File(eomip_h5name, 'w') as fout:
    fout.create_dataset('eip', data=eip)
    fout.create_dataset('vip', data=vip)
    '''.format(nroots, base_name))

def write_eomea(openfile, nroots, base_name):
    openfile.write('''
##############################
# EOM-EA-KRCCSD
##############################

# number of roots requested
# index(indices) of targetted k point(s) 
kptlist = [0]

myeom = eom_kccsd_rhf.EOMEA(mycc)
eea, vea = myeom.eaccsd(nroots={}, imds=imds, kptlist=kptlist)

eomea_h5name = '{}' + '_eomea.h5'
eomea_h5name = os.path.join(prefix, eomea_h5name)

print('Saving eea, vea to {{}}'.format(eomea_h5name))
with h5py.File(eomea_h5name, 'w') as fout:
    fout.create_dataset('eea', data=eea)
    fout.create_dataset('vea', data=vea)
    '''.format(nroots, base_name))

def slurm_submit(cluster, run_dir, run_name, slurm_file):
    if cluster == burg:
        os.chdir(run_dir)
        p = subprocess.run(['sbatch', '-J', run_name, slurm_file], capture_output=True)
    elif cluster == moto:
        p = subprocess.run(['ssh', '-t', 'terremoto', 'cd {} ; sbatch -J {} {}'.format(run_dir, run_name, slurm_file)], capture_output=True)
    return p.stdout.decode('utf-8').split()[-1]

def slurm_submit_dependency(cluster, run_dir, run_name, slurm_file, slurm_output):
    if cluster == burg:
        os.chdir(run_dir)
        p = subprocess.run(['sbatch', '--dependency=afterok:{}'.format(slurm_output), '-J', run_name, slurm_file], capture_output=True)
    elif cluster == moto:
        p = subprocess.run(['ssh', '-t', 'terremoto', 'cd {} ; sbatch --dependency=afterok:{} -J {} {}'.format(run_dir, slurm_output, run_name, slurm_file)], capture_output=True)
    return p.stdout.decode('utf-8').split()[-1]

def full_run(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, scaled_center, vb_nroots, cb_nroots, write_ip=True, write_ea=True):
    with open(run_file, 'w') as rf:
        write_header(rf)
        write_cell(rf, formula, scaled_center, basis_file)
        write_scf(rf, kdensity, base_name)
        write_ccsd(rf, base_name)
        if write_ip:
            write_eomip(rf, vb_nroots, base_name)
        if write_ea:
            write_eomea(rf, cb_nroots, base_name)
    p = slurm_submit(cluster, run_dir, run_name, slurm_file)
    return p
        

def full_ccsd(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, scaled_center, slurm_output=None, make_imds_ip=True, make_imds_ea=True):
    with open(run_file, 'w') as rf:
        write_header(rf)
        write_cell(rf, formula, scaled_center, basis_file)
        write_scf(rf, kdensity, base_name)
        write_ccsd(rf, base_name)
        if make_imds_ip:
            write_imds_ip(rf)
        if make_imds_ea:
            write_imds_ea(rf)
        if make_imds_ip or make_imds_ea:
            save_imds(rf)
    if slurm_output == None:
        p = slurm_submit(cluster, run_dir, run_name, slurm_file)
    else:
        p = slurm_submit_dependency(cluster, run_dir, run_name, slurm_file, slurm_output)
    return p

def run_eom(cluster, run_dir, run_name, run_file, base_name, slurm_file, slurm_output, formula, basis_file, kdensity, scaled_center, vb_nroots, cb_nroots, made_imds=True, write_ip=True, write_ea=True):
    with open(run_file, 'w') as rf:
        write_header(rf)
        write_cell(rf, formula, scaled_center, basis_file)
        write_scf(rf, kdensity, base_name)
        read_ccsd(rf, base_name)
        if made_imds:
            read_imds(rf)
        if write_ip:
            write_eomip(rf, vb_nroots, base_name)
        if write_ea:
            write_eomea(rf, cb_nroots, base_name)
    p = slurm_submit_dependency(cluster, run_dir, run_name, slurm_file, slurm_output)
    return p

def make_imds(cluster, run_dir, run_name, run_file, base_name, slurm_file, slurm_output, formula, basis_file, kdensity, scaled_center, make_imds_ip=True, make_imds_ea=True):
    with open(run_file, 'w') as rf:
        write_header(rf)
        write_cell(rf, formula, scaled_center, basis_file)
        write_scf(rf, kdensity, base_name)
        read_ccsd(rf, base_name)
        if make_imds_ip:
            write_imds_ip(rf)
        if make_imds_ea:
            write_imds_ea(rf)
        if make_imds_ip or make_imds_ea:
            save_imds(rf)
    p = slurm_submit_dependency(cluster, run_dir, run_name, slurm_file, slurm_output)
    return p

def cleanup(cluster, run_dir, run_name, run_file, base_name, slurm_file, slurm_output):
    with open(run_file, 'w') as rf:
        rf.write('''#!/usr/bin/env python
import os

run_dir = '{}'
base_name = '{}'
cderi_filename = base_name + '_cderi.h5'
cderi_filename = os.path.join(run_dir, cderi_filename)
cc_filename = base_name + '_cc.h5'
cc_filename = os.path.join(run_dir, cc_filename)

if os.path.isfile(cderi_filename):
    os.remove(cderi_filename)

if os.path.isfile(cc_filename):
    os.remove(cc_filename)
                '''.format(run_dir, base_name))
    p = slurm_submit_dependency(cluster, run_dir, run_name, slurm_file, slurm_output)
    return p

materials = ['c', 'si', 'sic']

vb_scaled_centers = {
        'c':[0., 0., 0.], 
        'si':[0., 0., 0.], 
        'sic':[0., 0., 0.], 
        'bn':[0., 0., 0.], 
        'bp':[0., 0., 0.], 
        'aln':[0., 0., 0.], 
        'alp':[0., 0., 0.], 
        'mgo':[0., 0., 0.], 
        'mgs':[0., 0., 0.], 
        'lih':[0.5, 0.,  0.5], 
        'lif':[0., 0., 0.], 
        'licl':[0.0530303,  0.0530303,  0.10606061]}

cb_scaled_centers = {
        'c':[0.36363636, 0.,         0.36363636], 
        'si':[0.3989899, 0.,        0.3989899], 
        'sic':[0.5, 0.,  0.5], 
        'bn':[0.5, 0.,  0.5], 
        'bp':[0.38888889, 0.,         0.38888889], 
        'aln':[0.0199005, 0.0199005, 0.       ], 
        'alp':[0.5, 0.,  0.5], 
        'mgo':[0., 0., 0.], 
        'mgs':[0.5, 0.,  0.5], 
        'lih':[0.5, 0.,  0.5], 
        'lif':[0., 0., 0.], 
        'licl':[0., 0., 0.]}

vb_nroots = {
        'c':3, 
        'si':3, 
        'sic':3, 
        'bn':3, 
        'bp':3, 
        'aln':1, 
        'alp':3, 
        'mgo':3, 
        'mgs':3, 
        'lih':1, 
        'lif':3, 
        'licl':1}

cb_nroots = {
        'c':1, 
        'si':1, 
        'sic':1, 
        'bn':1, 
        'bp':1, 
        'aln':1, 
        'alp':1, 
        'mgo':1, 
        'mgs':1, 
        'lih':1, 
        'lif':1, 
        'licl':1}

basis_sets = ['cc-pvdz-lc.dat', 'cc-pvtz-lc.dat', 'cc-pvqz-lc.dat']
basis_set_names = ['gth-cc-pvdz', 'gth-cc-pvtz', 'gth-cc-pvqz']

burg_set = set(['c', 'si', 'bn', 'bp', 'alp', 'mgo', 'lih', 'lif'])
moto_set = set(['sic', 'aln', 'mgs', 'licl'])

for basis_index, basis in enumerate(basis_sets):

    for kdensity in range(1, 5):

        for formula_index, formula in enumerate(materials):

            if formula in burg_set:
                cluster = burg
            else:
                cluster = moto

            basis_file = os.path.join(cluster, 'builds/ccgto/basis/gth-hf-rev', basis)
            slurm_dir = os.path.join(cluster, 'Scripts/pyscf')
            slurm_file = 'job-pyscf-highmem-node-long-SLURM.sh'
            slurm_file = os.path.join(slurm_dir, slurm_file)

            if vb_scaled_centers[formula] == cb_scaled_centers[formula]:

                base_name = formula + '_' + basis_set_names[basis_index] + '_' + str(kdensity) + str(kdensity) + str(kdensity)
                run_dir = os.path.join(cluster, 'Documents/benchmarking_ccgto', base_name)
                if not os.path.exists(run_dir):
                    os.makedirs(run_dir)
                
                if kdensity < 4:
                    run_name = base_name + '_ccsd'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p1 = full_ccsd(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, vb_scaled_centers[formula], slurm_output=None, make_imds_ip=False, make_imds_ea=False)

                    run_name = base_name + '_imds'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p2 = make_imds(cluster, run_dir, run_name, run_file, base_name, slurm_file, p1, formula, basis_file, kdensity, vb_scaled_centers[formula], make_imds_ip=True, make_imds_ea=True)

                    run_name = base_name + '_eomip'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p3 = run_eom(cluster, run_dir, run_name, run_file, base_name, slurm_file, p2, formula, basis_file, kdensity, vb_scaled_centers[formula], vb_nroots[formula], cb_nroots[formula], made_imds=True, write_ip=True, write_ea=False)

                    run_name = base_name + '_eomea'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p4 = run_eom(cluster, run_dir, run_name, run_file, base_name, slurm_file, p3, formula, basis_file, kdensity, cb_scaled_centers[formula], vb_nroots[formula], cb_nroots[formula], made_imds=True, write_ip=False, write_ea=True)

                    run_name = base_name + '_cleanup'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p5 = cleanup(cluster, run_dir, run_name, run_file, base_name, slurm_file, p4)

                else:
                    run_name = base_name + '_ccsd1'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p1 = full_ccsd(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, vb_scaled_centers[formula], make_imds_ip=False, make_imds_ea=False)

                    run_name = base_name + '_ccsd2'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p2 = full_ccsd(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, vb_scaled_centers[formula], slurm_output=p1, make_imds_ip=False, make_imds_ea=False)

                    run_name = base_name + '_imds'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p3 = make_imds(cluster, run_dir, run_name, run_file, base_name, slurm_file, p2, formula, basis_file, kdensity, vb_scaled_centers[formula], make_imds_ip=True, make_imds_ea=True)

                    run_name = base_name + '_eomip'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p4 = run_eom(cluster, run_dir, run_name, run_file, base_name, slurm_file, p3, formula, basis_file, kdensity, vb_scaled_centers[formula], vb_nroots[formula], cb_nroots[formula], made_imds=True, write_ip=True, write_ea=False)

                    run_name = base_name + '_eomea'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p5 = run_eom(cluster, run_dir, run_name, run_file, base_name, slurm_file, p4, formula, basis_file, kdensity, cb_scaled_centers[formula], vb_nroots[formula], cb_nroots[formula], made_imds=True, write_ip=False, write_ea=True)

                    run_name = base_name + '_cleanup'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p6 = cleanup(cluster, run_dir, run_name, run_file, base_name, slurm_file, p5)

            else:
                # Valence Band
                base_name = formula + '_' + basis_set_names[basis_index] + '_' + str(kdensity) + str(kdensity) + str(kdensity) + '_vb'
                run_dir = os.path.join(cluster, 'Documents/benchmarking_ccgto', base_name)
                if not os.path.exists(run_dir):
                    os.makedirs(run_dir)
                
                if kdensity < 4:
                    run_name = base_name + '_ccsd'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p1 = full_ccsd(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, vb_scaled_centers[formula], slurm_output=None, make_imds_ip=False, make_imds_ea=False)

                    run_name = base_name + '_imds'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p2 = make_imds(cluster, run_dir, run_name, run_file, base_name, slurm_file, p1, formula, basis_file, kdensity, vb_scaled_centers[formula], make_imds_ip=True, make_imds_ea=False)

                    run_name = base_name + '_eomip'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p3 = run_eom(cluster, run_dir, run_name, run_file, base_name, slurm_file, p2, formula, basis_file, kdensity, vb_scaled_centers[formula], vb_nroots[formula], cb_nroots[formula], made_imds=True, write_ip=True, write_ea=False)

                    run_name = base_name + '_cleanup'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p4 = cleanup(cluster, run_dir, run_name, run_file, base_name, slurm_file, p3)

                else:
                    run_name = base_name + '_ccsd1'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p1 = full_ccsd(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, vb_scaled_centers[formula], slurm_output=None, make_imds_ip=False, make_imds_ea=False)

                    run_name = base_name + '_ccsd2'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p2 = full_ccsd(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, vb_scaled_centers[formula], slurm_output=p1, make_imds_ip=False, make_imds_ea=False)

                    run_name = base_name + '_imds'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p3 = make_imds(cluster, run_dir, run_name, run_file, base_name, slurm_file, p2, formula, basis_file, kdensity, vb_scaled_centers[formula], make_imds_ip=True, make_imds_ea=False)

                    run_name = base_name + '_eomip'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p4 = run_eom(cluster, run_dir, run_name, run_file, base_name, slurm_file, p3, formula, basis_file, kdensity, vb_scaled_centers[formula], vb_nroots[formula], cb_nroots[formula], made_imds=True, write_ip=True, write_ea=False)

                    run_name = base_name + '_cleanup'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p5 = cleanup(cluster, run_dir, run_name, run_file, base_name, slurm_file, p4)

                # Conduction Band
                base_name = formula + '_' + basis_set_names[basis_index] + '_' + str(kdensity) + str(kdensity) + str(kdensity) + '_cb'
                run_dir = os.path.join(cluster, 'Documents/benchmarking_ccgto', base_name)
                if not os.path.exists(run_dir):
                    os.makedirs(run_dir)
                
                if kdensity < 4:
                    run_name = base_name + '_ccsd'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p1 = full_ccsd(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, cb_scaled_centers[formula], slurm_output=None, make_imds_ip=False, make_imds_ea=False)

                    run_name = base_name + '_imds'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p2 = make_imds(cluster, run_dir, run_name, run_file, base_name, slurm_file, p1, formula, basis_file, kdensity, cb_scaled_centers[formula], make_imds_ip=False, make_imds_ea=True)

                    run_name = base_name + '_eomea'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p3 = run_eom(cluster, run_dir, run_name, run_file, base_name, slurm_file, p2, formula, basis_file, kdensity, cb_scaled_centers[formula], vb_nroots[formula], cb_nroots[formula], made_imds=True, write_ip=False, write_ea=True)

                    run_name = base_name + '_cleanup'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p4 = cleanup(cluster, run_dir, run_name, run_file, base_name, slurm_file, p3)

                else:
                    run_name = base_name + '_ccsd1'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p1 = full_ccsd(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, cb_scaled_centers[formula], slurm_output=None, make_imds_ip=False, make_imds_ea=False)

                    run_name = base_name + '_ccsd2'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p2 = full_ccsd(cluster, run_dir, run_name, run_file, base_name, slurm_file, formula, basis_file, kdensity, cb_scaled_centers[formula], slurm_output=p1, make_imds_ip=False, make_imds_ea=False)

                    run_name = base_name + '_imds'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p3 = make_imds(cluster, run_dir, run_name, run_file, base_name, slurm_file, p2, formula, basis_file, kdensity, cb_scaled_centers[formula], make_imds_ip=False, make_imds_ea=True)

                    run_name = base_name + '_eomea'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p4 = run_eom(cluster, run_dir, run_name, run_file, base_name, slurm_file, p3, formula, basis_file, kdensity, cb_scaled_centers[formula], vb_nroots[formula], cb_nroots[formula], made_imds=True, write_ip=False, write_ea=True)

                    run_name = base_name + '_cleanup'
                    run_file = run_name + '.py'
                    run_file = os.path.join(run_dir, run_file)
                    p5 = cleanup(cluster, run_dir, run_name, run_file, base_name, slurm_file, p4)
