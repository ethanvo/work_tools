''' Driver for a simple CISD calculation.'''
import json
import h5py 
from numpy import ones
from pyscf.mcscf import CASSCF
from pyscf.scf import RHF
from pyscf.lib import chkfile as chk
from pyscf.mcscf.avas import avas

def main():
  '''Command-line usage.'''

  with open("mccalc.in.json",'r') as inpf:
    args = json.load(inpf)

  chkfile = args.pop('chkfile')
  mol = chk.load_mol(chkfile)
  mf = RHF(mol)
  mf.__dict__.update(chk.load(chkfile,'scf'))

  mccalc(mol, mf, **args)

def export_mc(mf, mc, mctype, energy, gdf=None, extra_data=None, to_afqmc=False):
  ''' Writing MC information to disk.'''
  if extra_data is None: extra_data = {}
  data = {
      'e_tot': energy, 
      'chkfile': mc.chkfile,
      'mo_coeff': mc.mo_coeff,
      'mctype': mctype,
      **extra_data
    }

  with h5py.File("mcdata.h5",'w') as outf:
    for key in data:
      outf[key] = data[key]

def mccalc(mol, mf, mctype, mcargs, nroots=1):
  ''' Interface between workflow and PySCF MC calculations.'''
  assert mctype.lower()=='casscf', "Only CASSCF maintained now."
  print(f"\nStarting {mctype} calculation.")
  casscfcalc(mol, mf, mcargs, nroots)
  print("\nDone.")

def casscfcalc(mol, mf, mcargs, nroots):
  state = mcargs.pop('state')

  if 'use_avas' in mcargs and mcargs.pop('use_avas'):
    ncas, nelecas, orbs = avas(mf, ['V 3d'], threshold=mcargs.pop('avas_thresh'), canonicalize=False)
  else:
    ncas = mcargs.pop('ncas')
    nelecas = mcargs.pop('nelecas')
    orbs = None

  if state == 'average':
    energy, ci = [], []
    mc = CASSCF(mf, nelecas=nelecas, ncas=ncas, **mcargs).state_average_(ones(nroots)/nroots)
    mc.fcisolver.nroots = nroots
    mc.chkfile = "mccalc.chk"
    #mc.fix_spin_(ss=ss)
    mc.kernel(orbs)
    energy  += list(mc.e_states)
    ci      += mc.ci
    #sss     += [ss for i in range(nroots)]

    print("\nOutputting data.")
    export_mc(mf, mc, 'casscf', energy, extra_data=dict(ci=ci, nroots=nroots, nelecas=mc.nelecas, ncas=mc.ncas))

  elif state == 'specific':
    mcs = [CASSCF(mf, nelecas=nelecas, ncas=ncas, **mcargs).state_specific_(i) for i in range(nroots)]
    energy, cis = [], []
    for idx, mc in enumerate(mcs):
      print(f"\nWorking on state specific = {idx}")
      mc.max_cycle_macro = 100
      mc.kernel(orbs)
      energy.append(mc.e_tot)
      cis.append(mc.ci)

    print("\nOutputting data.")
    export_mc(mf, mc, 'casscf', energy, extra_data=dict(ci=cis, nroots=nroots, nelecas=mc.nelecas, ncas=mc.ncas))
  
  return mc

if __name__=='__main__':
  main()
