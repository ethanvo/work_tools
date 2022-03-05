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
    au_to_ang = 0.5291772083

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

    extent = 4.0
    xmin = x.min() * ang_to_au - extent
    xmax = x.max() * ang_to_au + extent
    ymin = y.min() * ang_to_au - extent
    ymax = y.max() * ang_to_au + extent
    zmin = z.min() * ang_to_au - extent
    zmax = z.max() * ang_to_au + extent

    with open(basename + '_mep.inp', 'w') as fp:
        fp.write('{0:d}\n'.format(npoints**3))
        for ix in np.linspace(xmin, xmax, npoints, True):
            for iy in np.linspace(ymin, ymax, npoints, True):
                for iz in np.linspace(zmin, zmax, npoints, True):
                    fp.write('{0:12.6f} {1:12.6f} {2:12.6f}\n'.format(ix, iy, iz))

    subprocess.check_call(['orca_vpot', basename + '.gbw', basename + '.scfp',
            basename + '_mep.inp', basename + '_mep.out'])

    vpot = read_vpot(basename + '_mep.out')

    with open(basename + '_mep.cube', 'w') as fp:
        fp.write('Generated with ORCA\n')
        fp.write('Electrostatic potential for ' + basename + '\n')
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            len(atoms), xmin, ymin, zmin))
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            npoints, (xmax - xmin) / float(npoints - 1), 0.0, 0.0))
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            npoints, 0.0, (ymax - ymin) / float(npoints - 1), 0.0))
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            npoints, 0.0, 0.0, (zmax - zmin) / float(npoints - 1)))
        for i, atom in enumerate(atoms):
            index = elements.index(atom)
            fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}{4:12.6f}\n'.format(
                index, 0.0, x[i] * ang_to_au, y[i] * ang_to_au, z[i] * ang_to_au))

        m = 0
        n = 0
        vpot = np.reshape(vpot, (npoints, npoints, npoints))
        for ix in range(npoints):
            for iy in range(npoints):
                for iz in range(npoints):
                    fp.write('{0:14.5e}'.format(vpot[ix][iy][iz]))
                    m += 1
                    n += 1
                    if (n > 5):
                        fp.write('\n')
                        n = 0
                if n != 0:
                    fp.write('\n')
                    n = 0


    dx = (xmax - xmin) / float(npoints - 1)
    dy = (ymax - ymin) / float(npoints - 1)
    dz = (zmax - zmin) / float(npoints - 1)
    vpot = read_vpot(basename + '_mep.out')
    vpot = np.reshape(vpot, (npoints, npoints, npoints))
    efieldx, efieldy, efieldz = np.gradient(vpot, dx, dy, dz)
    normlist = np.zeros(1)
    efield = np.stack((efieldx, efieldy, efieldz), axis=-1)
    efieldnorm = np.linalg.norm(efield, axis=-1)

    with open(basename + '.out') as fp:
        for line in fp:
            if line.find("Total Dipole Moment    :") != -1:
                dipole_data = line.split()

    with open(basename + '_efield.cube', 'w') as fp:
        fp.write('Generated with ORCA\n')
        fp.write('Electrostatic potential for ' + basename + '\n')
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            len(atoms), xmin, ymin, zmin))
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            npoints, (xmax - xmin) / float(npoints - 1), 0.0, 0.0))
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            npoints, 0.0, (ymax - ymin) / float(npoints - 1), 0.0))
        fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(
            npoints, 0.0, 0.0, (zmax - zmin) / float(npoints - 1)))
        for i, atom in enumerate(atoms):
            index = elements.index(atom)
            fp.write('{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}{4:12.6f}\n'.format(
                index, 0.0, x[i] * ang_to_au, y[i] * ang_to_au, z[i] * ang_to_au))

        m = 0
        n = 0
        vpot = np.reshape(vpot, (npoints, npoints, npoints))
        for ix in range(npoints):
            for iy in range(npoints):
                for iz in range(npoints):
                    fp.write('{0:14.5e}'.format(efieldnorm[ix][iy][iz]))
                    m += 1
                    n += 1
                    if (n > 5):
                        fp.write('\n')
                        n = 0
                if n != 0:
                    fp.write('\n')
                    n = 0

    i = 0
    for ix in np.linspace(xmin, xmax, 10, endpoint=False):
        j = 0
        for iy in np.linspace(ymin, ymax, 10, endpoint=False):
            k = 0
            for iz in np.linspace(zmin, zmax, 10, endpoint=False):
                normlist = np.append(normlist, np.linalg.norm(np.array([efieldx[i,j,k],efieldy[i,j,k],efieldz[i,j,k]])))
                k += 10
            j += 10
        i += 10

    maxnorm = normlist.max()
    
    with open(basename + '-ESPwDipole_rotate.vmd', 'w') as fp:
        fp.write('mol new {' + basename + '.eldens.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 }\n')
        fp.write('mol addfile {' + basename + '_mep.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n')
        fp.write('mol addfile {' + basename + '_efield.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n')
        fp.write('color Display Background white\n')
        fp.write('mol modcolor 0 0 Element\n')
        fp.write('mol modstyle 0 0 Bonds 0.100000 12.000000\n')
        fp.write('mol color Element\n')
        fp.write('mol representation VDW 0.300000 12.000000\n')
        fp.write('mol selection all\n')
        fp.write('mol material Opaque\n')
        fp.write('mol addrep 0\n')
        fp.write('mol color Volume 1\n')
        fp.write('mol representation Isosurface 0.010000 0 0 0 1 1\n')
        fp.write('mol selection all\n')
        fp.write('mol material Transparent\n')
        fp.write('mol addrep 0\n')
        fp.write('mol scaleminmax 0 2' + '{0:12.6f}'.format(np.min(vpot)) + '{0:12.6f}'.format(np.max(vpot) * 0.0015) + '\n')
        fp.write('axes location Off\n')
        fp.write('scale by 0.8\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/vmd_draw_vector.tcl\n')
        fp.write('draw color black\n')
        fp.write('vmd_draw_vector 0 [list' + '{0:12.6f} {1:12.6f} {2:12.6f}'.format((xmin + xmax) / 2 * au_to_ang, (ymin + ymax) / 2 * au_to_ang, zmin) + ']')
        fp.write(' [list' + '{0:12.6f} {1:12.6f} {2:12.6f}'.format(float(dipole_data[4]), float(dipole_data[5]), float(dipole_data[6])) + '] 1.0 100 0.1\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/rotation_movie_ESP.tcl\n')
        fp.write('make_rotation_movie_files\n')
        fp.write('exit')

    with open(basename + '-FieldMap_rotate.vmd', 'w') as fp:
        fp.write('mol new {' + basename + '.eldens.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 }\n')
        fp.write('mol addfile {' + basename + '_mep.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n')
        fp.write('mol addfile {' + basename + '_efield.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n')
        fp.write('color Display Background white\n')
        fp.write('mol modcolor 0 0 Element\n')
        fp.write('mol modstyle 0 0 Bonds 0.100000 12.000000\n')
        fp.write('mol color Element\n')
        fp.write('mol representation VDW 0.300000 12.000000\n')
        fp.write('mol selection all\n')
        fp.write('mol material Opaque\n')
        fp.write('mol addrep 0\n')
        fp.write('axes location Off\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/vmd_tricolor_scale.tcl\n')
        fp.write('tricolor_scale\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/color_scale_bar_new.tcl\n')
        fp.write('namespace import ::ColorBar::*\n')
        fp.write('color_scale_bar 0.5 0.05 0 1 0 50 5\n')
        fp.write('scale by 0.8\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/vmd_draw_vector.tcl\n')
        i = 0
        for ix in np.linspace(xmin, xmax, 10, endpoint=False):
            j = 0
            for iy in np.linspace(ymin, ymax, 10, endpoint=False):
                k = 0
                for iz in np.linspace(zmin, zmax, 10, endpoint=False):
                    mag = np.linalg.norm(np.array([efieldx[i,j,k],efieldy[i,j,k],efieldz[i,j,k]]))
                    fp.write('draw color ')
                    if mag > 0.1:
                        fp.write('1056')
                    else:
                        fp.write(str(int(mag/0.097234519 * 1023 + 33)))
                    fp.write('\n')
                    fp.write('vmd_draw_vector 0 [list' + '{0:12.6f} {1:12.6f} {2:12.6f}'.format(ix * au_to_ang, iy * au_to_ang, iz * au_to_ang) + ']')
                    fp.write(' [list' + '{0:12.6f} {1:12.6f} {2:12.6f}'.format(efieldx[i,j,k] / mag, efieldy[i,j,k] / mag, efieldz[i,j,k] / mag) + '] 0.5 100 0.05\n')
                    k += 10
                j += 10
            i += 10
        fp.write('source /moto/berkelbach/users/eav2136/Applications/rotation_movie_FieldMap.tcl\n')
        fp.write('make_rotation_movie_files\n')
        fp.write('exit')

    with open(basename + '-FieldScan.vmd', 'w') as fp:
        fp.write('mol new {' + basename + '.eldens.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 }\n')
        fp.write('mol addfile {' + basename + '_mep.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n')
        fp.write('mol addfile {' + basename + '_efield.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n')
        fp.write('color Display Background white\n')
        fp.write('mol modcolor 0 0 Element\n')
        fp.write('mol modstyle 0 0 Bonds 0.100000 12.000000\n')
        fp.write('mol color Element\n')
        fp.write('mol representation VDW 0.300000 12.000000\n')
        fp.write('mol selection all\n')
        fp.write('mol material Opaque\n')
        fp.write('mol addrep 0\n')
        fp.write('axes location Off\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/vmd_tricolor_scale.tcl\n')
        fp.write('tricolor_scale\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/color_scale_bar_new.tcl\n')
        fp.write('namespace import ::ColorBar::*\n')
        fp.write('color_scale_bar 0.5 0.05 0 1 0 50 5\n')
        fp.write('scale by 0.8\n')
        fp.write('rotate y by -45\n')
        fp.write('mol color Volume 2\n')
        fp.write('mol representation VolumeSlice 0.000000 2.000000 2.000000 2.000000\n')
        fp.write('mol selection all\n')
        fp.write('mol material Opaque\n')
        fp.write('mol addrep 0\n')
        fp.write('mol scaleminmax 0 2 0.000000 0.097234519\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/slice_movie_FieldScan.tcl\n')
        fp.write('make_slice_movie_files\n')
        fp.write('exit')

    with open(basename + '-FieldMap_FieldScan.vmd', 'w') as fp:
        fp.write('mol new {' + basename + '.eldens.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 }\n')
        fp.write('mol addfile {' + basename + '_mep.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n')
        fp.write('mol addfile {' + basename + '_efield.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n')
        fp.write('color Display Background white\n')
        fp.write('mol modcolor 0 0 Element\n')
        fp.write('mol modstyle 0 0 Bonds 0.100000 12.000000\n')
        fp.write('mol color Element\n')
        fp.write('mol representation VDW 0.300000 12.000000\n')
        fp.write('mol selection all\n')
        fp.write('mol material Opaque\n')
        fp.write('mol addrep 0\n')
        fp.write('mol color Volume 1\n')
        fp.write('mol representation Isosurface 0.010000 0 0 0 1 1\n')
        fp.write('mol selection all\n')
        fp.write('mol material Transparent\n')
        fp.write('mol addrep 0\n')
        fp.write('mol scaleminmax 0 2' + '{0:12.6f}'.format(np.min(vpot)) + '{0:12.6f}'.format(np.max(vpot) * 0.0015) + '\n')
        fp.write('axes location Off\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/vmd_tricolor_scale.tcl\n')
        fp.write('tricolor_scale\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/color_scale_bar_new.tcl\n')
        fp.write('namespace import ::ColorBar::*\n')
        fp.write('color_scale_bar 0.5 0.05 0 1 0 50 5\n')
        fp.write('scale by 0.8\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/vmd_draw_vector.tcl\n')
        i = 0
        for ix in np.linspace(xmin, xmax, 10, endpoint=False):
            j = 0
            for iy in np.linspace(ymin, ymax, 10, endpoint=False):
                k = 0
                for iz in np.linspace(zmin, zmax, 10, endpoint=False):
                    mag = np.linalg.norm(np.array([efieldx[i,j,k],efieldy[i,j,k],efieldz[i,j,k]]))
                    fp.write('draw color ')
                    if mag > 0.1:
                        fp.write('1056')
                    else:
                        fp.write(str(int(mag/0.097234519 * 1023 + 33)))
                    fp.write('\n')
                    fp.write('vmd_draw_vector 0 [list' + '{0:12.6f} {1:12.6f} {2:12.6f}'.format(ix * au_to_ang, iy * au_to_ang, iz * au_to_ang) + ']')
                    fp.write(' [list' + '{0:12.6f} {1:12.6f} {2:12.6f}'.format(efieldx[i,j,k] / mag, efieldy[i,j,k] / mag, efieldz[i,j,k] / mag) + '] 0.5 100 0.05\n')
                    k += 10
                j += 10
            i += 10
        fp.write('source /moto/berkelbach/users/eav2136/Applications/rotation_movie_FieldMap_FieldScan.tcl\n')
        fp.write('make_rotation_movie_files\n')
        fp.write('mol color Volume 2\n')
        fp.write('mol representation VolumeSlice 0.000000 2.000000 2.000000 2.000000\n')
        fp.write('mol selection all\n')
        fp.write('mol material Opaque\n')
        fp.write('mol addrep 0\n')
        fp.write('mol scaleminmax 0 3 0.000000 0.097234519\n')
        fp.write('source /moto/berkelbach/users/eav2136/Applications/slice_movie_FieldMap_FieldScan.tcl\n')
        fp.write('make_slice_movie_files\n')
        fp.write('exit')