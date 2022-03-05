#!/usr/bin/env python
import os
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

burg = '/burg/berkelbach/users/eav2136'
moto = '/moto/berkelbach/users/eav2136'

au2ev = 27.211386245988
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

experiment = {'c':5.48, 'si':1.11, 'sic':2.42, 'bn':6.4, 'bp':2.4, 'alp':2.5, 'mgo':7.8, 'mgs':4.45, 'lih':4.8, 'lif':13.6, 'licl':10.4}

Nk_dz = np.arange(1, 4)
Nk_tz = np.arange(1, 3)
Nk_qz = np.arange(1, 2)
Nk_dz = 1/Nk_dz
Nk_tz = 1/Nk_tz
Nk_qz = 1/Nk_qz

eip = {}
eea = {}
bandgap = {}
for formula in materials:
    eip[formula] = np.zeros((3, 3))
    eea[formula] = np.zeros((3, 3))
    bandgap[formula] = np.zeros((3, 3))

for basis_index, basis in enumerate(basis_sets):

    for kdensity in range(1, 4):

        for formula_index, formula in enumerate(materials):

            if formula in burg_set:
                cluster = burg
            else:
                cluster = moto

            if vb_scaled_centers[formula] == cb_scaled_centers[formula]:

                base_name = formula + '_' + basis_set_names[basis_index] + '_' + str(kdensity) + str(kdensity) + str(kdensity)
                run_dir = os.path.join(cluster, 'Documents/benchmarking_ccgto', base_name)

                ip_h5name = os.path.join(run_dir, base_name + '_eomip.h5')
                ea_h5name = os.path.join(run_dir, base_name + '_eomea.h5')

                if os.path.exists(ip_h5name):
                    with h5py.File(ip_h5name, 'r') as fin:
                        print(fin['eip'][:])
                        eip[formula][kdensity - 1, basis_index] = np.amax(np.array(fin['eip'][:][:]))

                if os.path.exists(ea_h5name):
                    with h5py.File(ea_h5name, 'r') as fin:
                        eea[formula][kdensity - 1, basis_index] = np.amax(np.array(fin['eea'][:][:]))

            else:
                # Valence Band
                base_name = formula + '_' + basis_set_names[basis_index] + '_' + str(kdensity) + str(kdensity) + str(kdensity) + '_vb'
                run_dir = os.path.join(cluster, 'Documents/benchmarking_ccgto', base_name)

                ip_h5name = os.path.join(run_dir, base_name + '_eomip.h5')
                if os.path.exists(ip_h5name):
                    with h5py.File(ip_h5name, 'r') as fin:
                        print(fin['eip'][:])
                        eip[formula][kdensity - 1, basis_index] = np.amax(np.array(fin['eip'][:][:]))
                
                # Conduction Band
                base_name = formula + '_' + basis_set_names[basis_index] + '_' + str(kdensity) + str(kdensity) + str(kdensity) + '_cb'
                run_dir = os.path.join(cluster, 'Documents/benchmarking_ccgto', base_name)

                ea_h5name = os.path.join(run_dir, base_name + '_eomea.h5')
                if os.path.exists(ea_h5name):
                    with h5py.File(ea_h5name, 'r') as fin:
                        eea[formula][kdensity - 1, basis_index] = np.amax(np.array(fin['eea'][:][:]))
print(eip)
print(eea)

for formula in materials:
    bandgap[formula] = (eea[formula] + eip[formula]) * au2ev

print(bandgap)
for formula in materials:
    fig, ax = plt.subplots(figsize=(6, 5), dpi=100)

    ax.plot(Nk_dz, bandgap[formula][0:3, 0], color='c', marker='o', label='dz')
    ax.plot(Nk_tz, bandgap[formula][0:2, 1], color='m', marker='o', label='tz')
    ax.plot(Nk_qz, bandgap[formula][0:1, 2], color='y', marker='o', label='qz')
    ax.plot(0, experiment[formula.lower()], color='r', marker='o')

    ax.set_ylabel('Bandgap (eV)')
    ax.set_xlim(0)
    #ax.set_xticks(x_coords)
    #ax.set_xticklabels([r'$1^{-3}$', r'$2^{-3}$', r'$3^{-3}$'])

    plt.xlabel(r'$N_k^{-\frac{1}{3}}$')
    plt.title(formula)
    plt.tight_layout()
    plt.legend()

    plt.savefig(formula.lower() + '_preliminary_nk.png')
print('Script Finished!')
