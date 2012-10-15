import os
import re

import Logger

def SummarizePlots(PlotNames=list,
                   PlotSummaryBaseName=str,
                   Extension='.pdf',
                   Log=Logger):
    LogString = '** Summarizing plots to \"'+PlotSummaryBaseName+'.tex\" => \"'+PlotSummaryBaseName+'.pdf\" ..'
    print LogString
    Log.Write(LogString+'\n')
    DataFileNameListReverseCopy = []
    DataFileNameListReverseCopy.extend(PlotNames)
    DataFileNameListReverseCopy.reverse()
    fw = open(PlotSummaryBaseName+'.tex','w')
    fw.write('\documentclass[pre,amsmath,onecolumn,floatfix,showpacs,fleqn,a4paper,superscriptaddress]{revtex4}\n')
    fw.write('\usepackage{latexsym}\n')
    fw.write('\usepackage{amsmath,amsfonts,amssymb}\n')
    fw.write('\usepackage[dvips]{graphicx}\n')
    fw.write('\usepackage{subfigure}\n')
    fw.write('\usepackage{float}\n')
    fw.write('\usepackage{pifont}\n')
    fw.write('\usepackage{color}\n')
    fw.write('\n')
    fw.write('\\renewcommand{\\textfraction}{0.00} \\renewcommand{\\topfraction}{1.0}\n')
    fw.write('\\renewcommand{\\bottomfraction}{1.0}\n')
    fw.write('\\renewcommand{\\floatpagefraction}{1}\n')
    fw.write('\\newcommand\\vcent[1]{\ensuremath{\\vcenter{\hbox{{#1}}}}}\n')
    fw.write('\n')
    fw.write('\\begin{document}\n')
    NColumns     = 2
    NRows        = 5
    CountRows = 0
    while(len(DataFileNameListReverseCopy)>0):
        if(CountRows==0):
            fw.write('\\begin{center}\n')
            fw.write('\\begin{tabular}{cc}\n')
        if(len(DataFileNameListReverseCopy)==1):
            fw.write('\includegraphics[width=8cm]{'+
                     re.sub('.dat',Extension,DataFileNameListReverseCopy.pop())+'} \\\\\n')
            CountRows += 1
        else:
            fw.write('\includegraphics[width=8cm]{'+
                     re.sub('.dat',Extension,DataFileNameListReverseCopy.pop())+'} &\n'+
                     '\includegraphics[width=8cm]{'+
                     re.sub('.dat',Extension,DataFileNameListReverseCopy.pop())+'} \\\\\n')
            CountRows += 1
        if(CountRows==NRows):
            fw.write('\end{tabular}\n')
            fw.write('\end{center}\n')
            CountRows = 0
        if(CountRows==0):
            fw.write('\\newpage\n')
    if(CountRows!=0):
        fw.write('\end{tabular}\n')
        fw.write('\end{center}\n')
    fw.write('\end{document}\n')
    fw.close()
    del fw

    os.system('pdflatex '+PlotSummaryBaseName+'.tex > TeX.out 2>&1')

    return