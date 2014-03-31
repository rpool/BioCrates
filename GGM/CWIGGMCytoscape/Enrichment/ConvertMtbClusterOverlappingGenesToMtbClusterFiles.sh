#! /bin/bash

#tail -n +2 MtbClusterOverlappingGenes.csv | awk -F, '{print $1,$2}'
NClusters=`wc -l MtbClusterOverlappingGenes.csv | awk '{print $1-1}'`
for (( i=1; i<=$NClusters; i++ ))
do
    Line=`awk -F, '{print $1","$2}' MtbClusterOverlappingGenes.csv | grep $i\,`
    echo Cluster$i.txt $Line
    echo $Line | awk -F, '{print $2}' | sed -e s/\;/\\n/g > Cluster$i.txt
done

exit
