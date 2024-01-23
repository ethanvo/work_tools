#!/usr/bin/env python
import json
import os

def load(filename):
    with open(filename, 'r') as fin:
        data = json.loads(fin.read())
    return data

def dump(data, filename):
    with open(filename, 'w') as fout:
        fout.write(json.dumps(data, indent=4))

def remove(filename):
    if os.path.exists(filename):
        os.remove(filename)
    else:
        print("The file does not exist")
