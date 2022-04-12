import ase
import json
import h5py 
from os.path import exists
from os import rename
import pyscf
from numpy import eye
from pyscf.scf import fast_newton
from pyscf.gto import Mole
from pyscf.pbc.gto import Cell
from pyscf.pbc.tools.pyscf_ase import ase_atoms_to_pyscf

def main():
  ''' Command line interface.'''

  with open("scfcalc.in.json",'r') as inpf:
    args = json.load(inpf)

  scfcalc(**args)
  print("\nDone with scfcalc.")

def scfcalc(struct_source, structargs, mfargs, stability_check=False):
  '''Perform SCF for a given structure. '''
  print("Entering scfcalc.")

  mftype = mfargs.pop('mftype')

  struct = makestruct(struct_source, **structargs)
  mfcall = {
      ('RHF'  , False) : pyscf.scf.RHF,
      ('RKS'  , False) : pyscf.scf.RKS,
      ('ROHF' , False) : pyscf.scf.ROHF,
      ('ROKS' , False) : pyscf.scf.ROKS,
      ('UHF'  , False) : pyscf.scf.UHF,
      ('UKS'  , False) : pyscf.scf.UKS,
      ('RHF'  , True)  : pyscf.pbc.scf.RHF,
      ('RKS'  , True)  : pyscf.pbc.scf.RKS,
      ('ROHF' , True)  : pyscf.pbc.scf.ROHF,
      ('ROKS' , True)  : pyscf.pbc.scf.ROKS,
      ('UHF'  , True)  : pyscf.pbc.scf.UHF,
      ('UKS'  , True)  : pyscf.pbc.scf.UKS,
    }[mftype, 'a' in struct.__dict__]

  mf = mfcall(struct)
  mf.chkfile = "running.chk"

  if 'use_newton' in mfargs and mfargs.pop('use_newton'):
    print("Using Newton solver.")
    mf = fast_newton(mf)
  if 'usedf' in mfargs and mfargs.pop('usedf'):
    print("Using GDF.")
    mf = mf.density_fit()
    mf.with_df._cderi_to_save = "scfcalc_gdf.h5"

  if 'xc' in mfargs and mfargs['xc'].lower() != 'hf':
    mf.xc = mfargs.pop('xc')

  mf.max_cycle = 100
  mf.__dict__.update(**mfargs)
  mf.build()

  if exists("guess.chk"):
    print("Initiating from 'guess.chk'.")
    dm = mf.from_chk("guess.chk")
    mf.kernel(dm)
  else:
    mf.kernel()

  with h5py.File(mf.chkfile,'a') as outf:
    outf['scf/converged'] = mf.converged

  if mf.converged:
    print("Calculation converged, saving result to 'scfcalc.chk'.")
    rename(mf.chkfile, "scfcalc.chk")

  if stability_check:
    mf.stability()

ecp_basis_format = {
   'ccecp': 'ccecpccpv{basis}z',
   'bfdec': 'bfdv{basis}z',
   'none': 'ccpv{basis}z',
   'nonedkh': 'ccpv{basis}zdkh',
 }

def makestruct(struct_source, charge, multiplicity, basis, ecp, latc=None):
  ''' PySCF Mole object construction + apply some computation parameters. '''
  struct_ase = ase.io.read(struct_source)

  if latc is None:
    print("OBC construction.")
    struct = Mole()
  else:
    print("PBC construction.")
    struct = Cell()
    struct.a = latc * eye(3)

  struct.atom = ase_atoms_to_pyscf(struct_ase)
  struct.max_memory = 760e3
  struct.dimension = 3
  struct.symmetry = True

  struct.charge = charge
  struct.spin = multiplicity - 1
  struct.verbose = 4
  struct.ecp = ecp if 'none' not in ecp else None
  struct.basis = ecp_basis_format[ecp].format(basis=basis)
  struct.build()

  return struct

if __name__=='__main__':
  main()
