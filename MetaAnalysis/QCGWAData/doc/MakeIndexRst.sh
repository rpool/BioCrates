#! /bin/bash

PATHAPPEND=`echo $PWD | sed -e s/"\/doc"//`
export PYTHONPATH=$PATHAPPEND:$PYTHONPATH

generate_modules.py -n QCGWAData -d ./ -f -t ~/workspace/BioCrates/MetaAnalysis/QCGWAData

cat header.rst          >  index.rst
echo ""                 >> index.rst
#cat modules.rst         >> index.rst
#echo ""                 >> index.rst
cat Invocation.rst      >> index.rst
echo ""                 >> index.rst
cat QCGWAData.txt       >> index.rst
echo ""                 >> index.rst
cat ArgumentParser.txt  >> index.rst
echo ""                 >> index.rst
cat Logger.txt          >> index.rst
echo ""                 >> index.rst
cat DataContainer.txt   >> index.rst
echo ""                 >> index.rst
#cat File.txt            >> index.rst
#echo ""                 >> index.rst
#cat Format.txt          >> index.rst
#echo ""                 >> index.rst
#cat Checks.txt          >> index.rst
#echo ""                 >> index.rst
#cat Merge.txt           >> index.rst
#echo ""                 >> index.rst
#cat HapMap.txt          >> index.rst
#echo ""                 >> index.rst
#cat Filters.txt         >> index.rst
#echo ""                 >> index.rst
#cat FilterFunction.txt  >> index.rst
#echo ""                 >> index.rst
#cat Plotting.txt        >> index.rst

make html

exit
