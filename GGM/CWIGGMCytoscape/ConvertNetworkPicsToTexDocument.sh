#! /bin/bash

if [ -e Body.tex ]; then
    rm Body.tex
fi

#inkscape -f legend.svg -E legend.eps
#epstopdf legend.eps

#for FILENAME in KORAGGMLegendPartialCorrelations.svg KORAGGMLegendJaccardIndices.svg KORAGGMLegendPearsonCorrelations.svg
#do
#    BASENAME=$(basename $FILENAME | awk -F. '{print $1}')
#    sed -e s/font-size:18/font-size:20/ \
#        -e s/font-size:14/font-size:18/ \
#        -e s/font-size:16/font-size:20/ \
#        -e s/font-size:12/font-size:18/ ${BASENAME}.svg > tmp.svg
#    inkscape -f tmp.svg -E ${BASENAME}.eps
#    epstopdf ${BASENAME}.eps
#    rm tmp.svg
#done

for TYPE in PartialCorrelations
do
    for i in {0..9}
    do
        Cluster=""
        if [ $i == 0 ]; then
            Cluster=""
        else
            Cluster=Cluster"$i"
        fi
        echo ${TYPE} $Cluster
        CLUSTERFNAME=KORAGGM${TYPE}${Cluster}.pdf
        CLUSTERFNAMEC=KORAGGM${TYPE}${Cluster}Cropped.pdf
        FLEGEND=KORAGGMLegend${TYPE}.pdf
        FCAPTION=KORAGGM${TYPE}${Cluster}.tex
#        pdfcrop ${CLUSTERFNAME} ${CLUSTERFNAMEC}
        echo "\\clearpage"                                                                                                 >> Body.tex
        echo "\\thispagestyle{empty}"                                                                                      >> Body.tex
        echo "\\vspace*{\\stretch{1}}"                                                                                     >> Body.tex
        echo "\\begin{center}"                                                                                             >> Body.tex
        echo "\\begin{tabular}{c|c}"                                                                                       >> Body.tex
        echo "\\vcent{\\includegraphics[height=14cm]{$CLUSTERFNAMEC}} & \\vcent{\\includegraphics[height=18cm]{$FLEGEND}}" >> Body.tex
        echo "\\end{tabular}"                                                                                              >> Body.tex
        echo "\\end{center}"                                                                                               >> Body.tex
        echo ""                                                                                                            >> Body.tex
        echo "\\vspace{0.5cm}"                                                                                             >> Body.tex
        cat ${FCAPTION}                                                                                                    >> Body.tex
        echo "\\vspace*{\\stretch{3}}"                                                                                     >> Body.tex
        echo "\\newpage"                                                                                                   >> Body.tex
        echo ""                                                                                                            >> Body.tex

    done
done

cat Header.tex  > SimilarityGroupsCWIOnGGM.tex
cat Body.tex   >> SimilarityGroupsCWIOnGGM.tex
cat Footer.tex >> SimilarityGroupsCWIOnGGM.tex

pdflatex SimilarityGroupsCWIOnGGM.tex
pdflatex SimilarityGroupsCWIOnGGM.tex

if [ -e Body.tex ]; then
    rm Body.tex
fi

#inkscape -f legend.svg -E legend.eps
#epstopdf legend.eps

#for FILENAME in KORAGGMLegendPartialCorrelations.svg KORAGGMLegendJaccardIndices.svg KORAGGMLegendPearsonCorrelations.svg
#do
#    BASENAME=$(basename $FILENAME | awk -F. '{print $1}')
#    sed -e s/font-size:18/font-size:20/ \
#        -e s/font-size:14/font-size:18/ \
#        -e s/font-size:16/font-size:20/ \
#        -e s/font-size:12/font-size:18/ ${BASENAME}.svg > tmp.svg
#    inkscape -f tmp.svg -E ${BASENAME}.eps
#    epstopdf ${BASENAME}.eps
#    rm tmp.svg
#done

for TYPE in PartialCorrelations JaccardIndices PearsonCorrelations
do
    for i in {0..9}
    do
        Cluster=""
        if [ $i == 0 ]; then
            Cluster=""
        else
            Cluster=Cluster"$i"
        fi
        echo ${TYPE} $Cluster
        CLUSTERFNAME=KORAGGM${TYPE}${Cluster}.pdf
        CLUSTERFNAMEC=KORAGGM${TYPE}${Cluster}Cropped.pdf
        FLEGEND=KORAGGMLegend${TYPE}.pdf
        FCAPTION=KORAGGM${TYPE}${Cluster}.tex
#        pdfcrop ${CLUSTERFNAME} ${CLUSTERFNAMEC}
        echo "\\clearpage"                                                                                                 >> Body.tex
        echo "\\thispagestyle{empty}"                                                                                      >> Body.tex
        echo "\\vspace*{\\stretch{1}}"                                                                                     >> Body.tex
        echo "\\begin{center}"                                                                                             >> Body.tex
        echo "\\begin{tabular}{c|c}"                                                                                       >> Body.tex
        echo "\\vcent{\\includegraphics[height=14cm]{$CLUSTERFNAMEC}} & \\vcent{\\includegraphics[height=18cm]{$FLEGEND}}" >> Body.tex
        echo "\\end{tabular}"                                                                                              >> Body.tex
        echo "\\end{center}"                                                                                               >> Body.tex
        echo ""                                                                                                            >> Body.tex
        echo "\\vspace{0.5cm}"                                                                                             >> Body.tex
        cat ${FCAPTION}                                                                                                    >> Body.tex
        echo "\\vspace*{\\stretch{3}}"                                                                                     >> Body.tex
        echo "\\newpage"                                                                                                   >> Body.tex
        echo ""                                                                                                            >> Body.tex

    done
done

cat Header.tex  > Document.tex
cat Body.tex   >> Document.tex
cat Footer.tex >> Document.tex

pdflatex Document.tex
pdflatex Document.tex

exit
