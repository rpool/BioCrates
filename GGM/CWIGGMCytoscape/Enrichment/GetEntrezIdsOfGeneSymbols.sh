#! /bin/bash

for i in Cluster?.txt
do
    echo $i
    FContent=`cat $i`
    for Line in $FContent
    do
        Grep=`grep $Line$'\t' UniqGeneSymbols.txt_1386689090.txt_complete_info_results.txt | awk -F"\\t" '{print $2}'`
        echo $Grep
    done > Entrez$i
done

exit
