
"""
(c) 2013, 2019 Marius Retegan
License: BSD-2-Clause
Description: Create a .cube file of the moleclar electrostatic
             potential (MEP) using ORCA.
Run: python mep.py basename npoints (e.g. python mep.py water 40).
Arguments: basename - file name without the extension;
                      this should be the same for the .gbw and .scfp.
           npoints  - number of grid points per side
                      (80 should be fine)
"""

#!/usr/bin/env python


def read_xyz(xyz):
    atoms = list()
    x = list()
    y = list()
    z = list()

    with open(xyz) as fp:
        # Skip the first two lines.
        next(fp)
        next(fp)
        for line in fp:
            data = line.split()
            atoms.append(data[0])
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))

    return atoms, np.array(x), np.array(y), np.array(z)


def read_vpot(vpot):
    v = list()

    with open(vpot) as fp:
        next(fp)
        for line in fp:
            data = line.split()
            v.append(float(data[3]))

    return np.array(v)

if __name__ == '__main__':
    import os
    import sys
    import subprocess
    import numpy as np

    ang_to_au = 1.0 / 0.5291772083

    elements = [None,
         'H', 'He',
         'Li', 'Be',
         'B', 'C', 'N', 'O', 'F', 'Ne',
         'Na', 'Mg',
         'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
         'K', 'Ca',
         'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
         'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
         'Rb', 'Sr',
         'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
         'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
         'Cs', 'Ba',
         'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
         'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
         'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
         'Fr', 'Ra',
         'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
         'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Uub']

    basename = sys.argv[1]
    xyz = basename + '.xyz'

    if not os.path.isfile(xyz):
        sys.exit('Could not find the .xyz. To quickly generate one for '
                 'your molecule run: echo 11 | orca_plot {}.gbw -i.'.format(basename))

    atoms, x, y, z = read_xyz(xyz)

    try:
        npoints = int(sys.argv[2])
    except ValueError:
        sys.exit('Invalid number of points: {}'.format(sys.argv[2]))

    natoms = len(atoms)

    extent = 1.0
    xmin = x.min() * ang_to_au - extent
    xmax = x.max() * ang_to_au + extent
    ymin = y.min() * ang_to_au - extent
    ymax = y.max() * ang_to_au + extent
    zmin = z.min() * ang_to_au - extent
    zmax = z.max() * ang_to_au + extent

    with open(basename + '_efield.inp', 'w') as fp:
        fp.write('{0:d}\n'.format(npoints**3))
        for ix in np.linspace(xmin, xmax, npoints, True):
            for iy in np.linspace(ymin, ymax, npoints, True):
                for iz in np.linspace(zmin, zmax, npoints, True):
                    fp.write('{0:12.6f} {1:12.6f} {2:12.6f}\n'.format(ix, iy, iz))

    subprocess.check_call(['orca_vpot', basename + '.gbw', basename + '.scfp',
            basename + '_efield.inp', basename + '_efield.out'])

    dx = (xmax - xmin) / float(npoints - 1)
    dy = (ymax - ymin) / float(npoints - 1)
    dz = (zmax - zmin) / float(npoints - 1)
    vpot = read_vpot(basename + '_efield.out')
    vpot = np.reshape(vpot, (npoints, npoints, npoints))
    efieldx, efieldy, efieldz = np.gradient(vpot, dx, dy, dz)
    efieldxlist = np.reshape(efieldx, npoints ** 3)
    efieldylist = np.reshape(efieldy, npoints ** 3)
    efieldzlist = np.reshape(efieldz, npoints ** 3)

    with open(basename + '_efield_triples.txt', 'w') as fp:
        i = 0
        for ix in np.linspace(xmin, xmax, 10, True):
            j = 0
            for iy in np.linspace(ymin, ymax, 10, True):
                k = 0
                for iz in np.linspace(zmin, zmax, 10, True):
                    fp.write('vmd_draw_vector 0 [list' + '{0:12.6f} {1:12.6f} {2:12.6f}'.format(ix, iy, iz) + ']')
		    mag = np.linalg.norm(np.array([efieldx[i,j,k],efieldy[i,j,k],efield[i,j,k]])) 
                    fp.write(' [list' + '{0:12.6f} {1:12.6f} {2:12.6f}'.format(efieldx[i,j,k] / mag, efieldy[i,j,k] / mag, efieldz[i,j,k] / mag) + '] 1 100 0.1\n')
                    k += 10
                j += 10
            i += 10
