#!/usr/bin/env python
import sys
import numpy as np

basename = sys.argv[1]
field = []
with open(basename + '.output') as fp:
    for line in fp:
        if line.find("The field strength is") != -1:
            line_data = line.split()
            field.append(float(line_data[4]))
field = np.array(field)
print('Average hydrogen field strength:', np.mean(field[0:4]), 'V/nm')
print('Average carbon field strength:', np.mean(field[4:]), 'V/nm')