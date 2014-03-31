#! /bin/bash

awk -F\\t '{print $2}' UniqGeneSymbols.txt_1386949612.txt_complete_info_results.txt | grep -v NULL > NULLFilteredUniqEntrezIds.txt

exit
