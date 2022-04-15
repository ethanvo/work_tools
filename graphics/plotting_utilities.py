#!/usr/bin/env python
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from cycler import cycler

# Set default sans-serif font to Arial
matplotlib.rcParams['font.sans-serif'] = 'Arial'

# Set default to sans-serif
matplotlib.rcParams['font.family'] = 'sans-serif'

# Set color palette
matplotlib.rcParams['axes.prop_cycle'] = cycler('color', ['#648FFF', '#DC267F', '#FFB000', '#785EF0', '#FE6100'])

# Set Font Size
matplotlib.rcParams['font.size'] = 16

# Set Line Width
matplotlib.rcParams['lines.linewidth'] = 2

# Set Axes Line Width
matplotlib.rcParams['axes.linewidth'] = 1

# Set Legend Frame Off
matplotlib.rcParams['legend.frameon'] = False

# Example plot
'''
fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=600)
ax.plot(r, zerofield_shift, label='CASSCF 0 V/nm', color='b')
ax.legend()
ax.set_xlabel(r'O-O Distance ($\AA$)')
ax.set_ylabel('Energy (eV)')

# Remove box around plot
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.tight_layout()

plt.savefig('dissociation_shift_dft.svg')
'''
