#!/bin/bash
##imax stands for total number of total number of structural elements+1: for instance 141 for a 140 amino acid long proten, 21 for a 20 base pair long DNA poly nucleotide etc
##in split command below 'split -l 100 --numeric-suffixes=1 --suffix-length=3 --additional-suffix="-*.dat"  *-smp.dat "" --verbose', the number 100 to be replaced by the number = nwindow*nconf

##in split command below 'split -l 50 --numeric-suffixes=1 --suffix-length=2 --additional-suffix="-apo-smp.dat"  ${i}-apo.dat "" --verbose', the number 50 to be replaced by the number = nconf

imax=141
split -l 100 --numeric-suffixes=1 --suffix-length=3 --additional-suffix="-apo.dat"  apo-smp.dat "" --verbose
split -l 100  --numeric-suffixes=1 --suffix-length=3 --additional-suffix="-holo.dat" holo-smp.dat "" --verbose
for (( i=1; i<10; i++ ))
do
    mv 00${i}-apo.dat  ${i}-apo.dat
    mv 00${i}-holo.dat ${i}-holo.dat
done
for (( i=10; i<100; i++ ))
do
    mv 0${i}-apo.dat  ${i}-apo.dat
    mv 0${i}-holo.dat ${i}-holo.dat
done
for (( i=1; i<imax; i++ ))
do
    mkdir ${i}
    mv ${i}-apo.dat ${i}-holo.dat ${i}
    cp listhistprotmaxmin  parameter_maxmin histogram_max_min.py max_min_merge.sh max_min.py num_window  listhistprot histogram.py histparameter conf-thermodynamics.f90 listconformational_td ${i}
done
for (( i=1; i<imax; i++ ))
do
    cd ${i}
    split -l 50 --numeric-suffixes=1 --suffix-length=2 --additional-suffix="-apo-smp.dat"  ${i}-apo.dat "" --verbose
    split -l 50 --numeric-suffixes=1 --suffix-length=2 --additional-suffix="-holo-smp.dat"  ${i}-holo.dat "" --verbose
    rm *apo.dat
    rm *holo.dat
    python3 histogram_max_min.py < parameter_maxmin
    chmod +x max_min_merge.sh
    ./max_min_merge.sh
    python3 max_min.py < num_window
    rm max_min_apo.dat
    rm max_min_holo.dat
    rm max_min.dat
    cat histparameter max_min_convert.dat > histparam
    python3 histogram.py < histparam
    gfortran conf-thermodynamics.f90
    ./a.out < listconformational_td
    mv conftherm_av.dat ../${i}-conf-therm_av.dat
    cd ../
done
for (( i=1; i<10; i++ ))
do
    mv ${i}-conf-therm_av.dat  00${i}-conf-therm_av.dat
done
for (( i=10; i<100; i++ ))
do
    mv ${i}-conf-therm_av.dat  0${i}-conf-therm_av.dat
done

    
cat *-conf-therm_av.dat > conf-final.dat 
rm *-conf-therm_av.dat

for (( i=1; i<imax; i++ ))
do
    rm -r ${i}

done
paste resname.dat conf-final.dat > conf-therm.dat

 rm listhistprotmaxmin  parameter_maxmin  num_window  listhistprot histparameter  listconformational_td conf-final.dat resname.dat
 rm apo-smp.dat holo-smp.dat
 rm a.out
exit;
