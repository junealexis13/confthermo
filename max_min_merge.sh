# For apo protein and holo protein, this script can be used to combine all maximum and minimum values of each window
#!/bin/bash
cat *apo-max_min.dat > max_min_apo.dat
cat *holo-max_min.dat > max_min_holo.dat
rm *apo-max_min.dat
rm *holo-max_min.dat
