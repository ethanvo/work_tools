#!/usr/bin/env python3
import numpy as np
from pyscf.pbc import gto, scf
from pyscf.pbc.tools import lattice, pyscf_ase
from ase.dft.kpoints import get_bandpath
from fileutils import dump
import h5py

au2ev = 27.211386245988

def make_cell(formula, basis, pseudo=None, exp_to_discard=None):
    cell = gto.Cell()
    ase_atom = lattice.get_ase_atom(formula)
    cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atom)
    cell.a = ase_atom.cell[:]
    cell.unit = "B"
    cell.basis = basis
    cell.pseudo = pseudo
    cell.exp_to_discard = exp_to_discard
    cell.verbose = 7
    cell.build()
    return cell

def get_all_electron_bandstructure(formula, cell, path, npoints, kmesh=[2, 2, 2], e_kn_file=None, output_file=None):
    ase_atom = lattice.get_ase_atom(formula)
    bandpath = get_bandpath(path, ase_atom.cell, npoints=npoints)
    band_kpts_scaled = bandpath.kpts
    e_kn = []
    dm = None
    for kcenter in band_kpts_scaled:
        kpts = cell.make_kpts(kmesh, scaled_center=kcenter)
        mymf = scf.KRHF(cell, kpts=kpts, exxdiv="ewald").density_fit()
        escf = mymf.kernel(dm0=dm)
        e_kn.append(mymf.mo_energy[0])
    vbmax = -np.inf
    vbmax_k = 0
    for idx, en in enumerate(e_kn):
        vb_k = en[cell.nelectron // 2 - 1]
        if vb_k > vbmax:
            vbmax = vb_k
            vbmax_k = idx
    e_kn = np.asarray([en - vbmax for en in e_kn])
    cbmin_k = np.asarray(np.unravel_index(np.argmin(np.where(e_kn > 0.0, e_kn, np.inf), axis=None), e_kn.shape))
    vbmax_k_array = e_kn[vbmax_k]
    cbmin_k_array = e_kn[cbmin_k[0]]
    g_vbmax = vbmax_k_array[np.abs(vbmax_k_array) < 0.001].size
    g_cbmin = cbmin_k_array[cbmin_k_array - np.amin(np.where(cbmin_k_array > 0.0, cbmin_k_array, np.inf)) < 0.001].size
    if e_kn_file is not None:
        with h5py.File(e_kn_file, "w") as f:
            f.create_dataset("e_kn", data=e_kn)
    if output_file is not None:
        dump({"vbmax_k": vbmax_k, "cbmin_k": cbmin_k, "g_vbmax": g_vbmax, "g_cbmin": g_cbmin}, output_file)
    vbmax_kpt = band_kpts_scaled[vbmax_k]
    cbmin_kpt = band_kpts_scaled[cbmin_k[0]]
    return vbmax_kpt, cbmin_kpt, g_vbmax, g_cbmin
