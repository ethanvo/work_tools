#!/bin/bash
salloc -t 5-00:00:00 -A berkelbach --exclusive --nodes=1 --ntasks-per-node=1 --cpus-per-task=32 --mem=180G --mail-type=ALL --mail-user=eav2136@columbia.edu
