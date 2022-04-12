#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt

# Set default sans-serif font to Arial
matplotlib.rcParams['font.sans-serif'] = 'Arial'

# Set default to sans-serif
matplotlib.rcParams['font.family'] = 'sans-serif'

# Accessible colors
color_1 = '#648FFF'
color_2 = '#DC267F'
color_3 = '#FFB000'
color_4 = '#785EF0'
color_5 = '#FE6100'

# Example plot

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
