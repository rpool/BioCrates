#! /bin/bash

awk '{print $1}' PINA_HSA-20121210.sif.txt   > tmp.txt
awk '{print $NF}' PINA_HSA-20121210.sif.txt >> tmp.txt

sort tmp.txt | uniq > UniqueUniprotIDsInPINA.txt
rm tmp.txt

exit
