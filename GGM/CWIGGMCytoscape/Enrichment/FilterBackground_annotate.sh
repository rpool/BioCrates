#! /bin/bash

CutLineNo=`grep -n "\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-" Background.annotate | head -1 | cut -d: -f1`

head -$CutLineNo Background.annotate | grep . | grep -v \# | grep -v None | awk '{print $1}' > NoneFilteresEntrezBackground_annotate.txt

exit
