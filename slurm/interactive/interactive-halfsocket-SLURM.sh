#!/bin/bash
salloc -t 0-12:00:00 -A berkelbach --nodes=1 --sockets-per-node=1 --cores-per-socket=8 --threads-per-core=1 --ntasks-per-node=1 --cpus-per-task=8 --mem=45G --mail-type=ALL --mail-user=eav2136@columbia.edu
