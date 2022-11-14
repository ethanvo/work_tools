#!/bin/bash
salloc -t 0-12:00:00 -A berkelbach --nodes=1 --sockets-per-node=1 --cores-per-socket=16 --threads-per-core=1 --ntasks-per-node=1 --cpus-per-task=16 --mem=90G --mail-type=ALL --mail-user=eav2136@columbia.edu
