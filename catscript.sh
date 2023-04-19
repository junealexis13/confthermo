
# This script can be used to make aseparate file 'all-freeenergychange.dat' and  'all-entropychange.dat'  
# Remove all unnecessary files
# file which contains all free energy and entropy change, obtained from various block
#!/bin/bash
cat *-free-encost-res.dat > all-freeenergychange.dat
cat *-entropy-encost-res.dat > all-entropycost.dat

rm -rf *-free-encost-res.dat
rm -rf *-entropy-encost-res.dat
rm -rf *apo-hist.dat
rm -rf *holo-hist.dat
rm -rf *apo.dat
rm -rf *holo.dat
