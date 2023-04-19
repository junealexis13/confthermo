# This script can be used to split the random generated data set of dihedral angle according to user's choice
# The good choice of data set is 12
# For 960 frame conformation one can made 12 data set, each set contains 80 values.
# change split command accordingly
#!/bin/bash
#split -l300 --numeric-suffixes=1 --suffix-length=1 --additional-suffix=".lst"  file ""

split -l 80 --numeric-suffixes=1 --suffix-length=2 --additional-suffix="-apo.dat"  apo-protein-random-structure.dat "" --verbose
split -l 80 --numeric-suffixes=1 --suffix-length=2 --additional-suffix="-holo.dat"  holo-protein-random-structure.dat "" --verbose

