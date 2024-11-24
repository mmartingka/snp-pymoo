#!/bin/bash

checkPIDs() {
    n=`echo $pid | wc -w`
    echo "$n"
    while [ "$n" -gt 7 ]; do
        pid2=""
        for p in $pid; do
          if ps -p $p > /dev/null
          then
            pid2="$pid2 $p"
          fi
        done
        pid=$pid2
        n=`echo $pid | wc -w`
        echo $n
    done
}


#Valores comunes para todas las ejecuciones
dir_db=/home/mery/Desktop/TFM/venv/tfm_v0/data
dir_res=/home/mery/Desktop/TFM/venv/tfm_v0/results/db_1000SNPs_4000pac
dir_main=/home/mery/Desktop/TFM/venv/tfm_v0/snp/main.py

lmut=10
fmut=30
probCross=50

# ********************************** MEDIUM DB **********************************
for dimEpi in 2 5 8; do
    for instancia in $(seq 0 11); do
        file=$dir_db/db_1000SNPs_4000pac.txt;
        nohup python3 $dir_main 100 1500 $lmut $fmut $file $dimEpi $probCross > $dir_res/epi_$dimEpi/SNP-UNSGA3-$instancia-$dimEpi-PF.csv &
        p="$!"
        pid="$pid $p"
        checkPIDs
    done
done
