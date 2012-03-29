#! /usr/bin/env python
import argparse
import os
import scipy
import pylab
import re

import Logger

def LogArguments(Log=Logger,
                 ArgParser=argparse.ArgumentParser,
                 Arguments=argparse.Namespace):
    ArgParser.print_help()
    ArgParser.print_help(Log.GetFileHandle())

#   Calculate max. of keylength for formatting
    MaxLen = 0
    for Key in vars(Arguments).iterkeys():
        MaxLen = max(MaxLen,len(Key))
    LogString =  '\n****\n'+\
                 'Used arguments:\n'+\
                 '---------------'
    print LogString
    Log.Write(LogString+'\n')

    FormatString = '{:<'+str(MaxLen)+'}'
    LogString    = ''
    for Key,Value in vars(Arguments).iteritems():
        LogString += FormatString.format(Key)+': '+str(Value)+'\n'
    LogString += '****\n'
    print LogString
    Log.Write(LogString+'\n')
    return

def ParseArguments(Log=None):
    ArgumentParser = argparse.ArgumentParser(description=\
                                             'This python module can be used for visualizing the correlations between metabolite pairs.')
    ArgumentParser.add_argument('-i',
                                '--ratiofile',
                                dest='RatioFileName',
                                help='FILENAME: Name of the file containing the correlation data between metabolite pairs (input)',
                                metavar='FILENAME')
    ArgumentParser.add_argument('-l',
                                '--logtransformed',
                                dest='boLogTransformed',
                                help='FLAG: Is the data in the correlation data file based on log transformed concentrations (ln{c+1.0})? (input)',
                                action='store_true',
                                default=False)

    Arguments = ArgumentParser.parse_args()

    return ArgumentParser,\
           Arguments

def PylabUpdateParams():
    figwidth_pt   = 246.0 # pt (from revtex \showthe\columnwidth)
    inches_per_pt = 1.0/72.27
    figwidth      = figwidth_pt*inches_per_pt
    golden_mean   = (scipy.sqrt(5.0)-1.0)/2.0 # Aesthetic ratio
    figheight     = figwidth*golden_mean
    fig_size      = [figwidth,figheight]
    params        = {'backend': 'pdf',
                     'patch.antialiased': True,
                     'axes.labelsize': 8,
                     'axes.linewidth': 0.5,
                     'grid.color': '0.75',
                     'grid.linewidth': 0.25,
                     'grid.linestyle': ':',
                     'axes.axisbelow': False,
                     'text.fontsize': 8,
                     'legend.fontsize': 5,
                     'xtick.labelsize': 8,
                     'ytick.labelsize': 8,
                     'text.usetex': True,
                     'figure.figsize': fig_size}
    left   = 0.16
    bottom = 0.16
    width  = 0.86-left
    height = 0.95-bottom

    return params,\
           [left, bottom, width, height]

def PlotHeatMap(X=scipy.array,
                Y=scipy.array,
                Z=scipy.array,
                XLabel='',
                YLabel='',
                ZLabel=''):

    PylabParameters,\
    Rectangle         = PylabUpdateParams()
    pylab.rcParams.update(PylabParameters)

    PylabFigure = pylab.figure()
    PylabFigure.clf()

    PylabAxis   = PylabFigure.add_axes(Rectangle)
    PylabAxis.set_xlabel(XLabel)
    PylabAxis.set_ylabel(YLabel)
    Values = pylab.ones(shape=(163,163),dtype=float) # indices run from 1 in data file!
    for z in range(len(Z)):
        i = X[z]
        j = Y[z]
        Values[i-1][j-1] = Z[z] # indices run from 1 in data file!
#    Extent      = [X[0],X[-1],Y[0],Y[-1]]
#    PylabAxis.matshow(Values,extent=Extent)
#    Extent        = [X[0]-.5,X[-1]-.5,Y[0]-.5,Y[-1]-.5]
    Extent        = [X.min()-0.5,X.max()+0.5,Y.min()-0.5,Y.max()+0.5]
    CMap          = pylab.cm.get_cmap(name='jet',
                                      lut=None)
    PylabImAxis   = PylabFigure.add_axes(Rectangle)
    PylabImage    = PylabImAxis.imshow(Values,
                                       extent=Extent,
                                       interpolation='nearest',
                                       cmap=CMap,
                                       origin='lower')

    PylabAxis.plot([41.0,41.0],
                   [1.0,41.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([1.0,1.0],
                   [1.0,41.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([1.0,41.0],
                   [41.0,41.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([1.0,41.0],
                   [1.0,1.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([55.0,55.0],
                   [41.0,55.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([41.0,41.0],
                   [41.0,55.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([41.0,55.0],
                   [55.0,55.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([41.0,55.0],
                   [41.0,41.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([147.0,147.0],
                   [55.0,147.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([55.0,55.0],
                   [55.0,147.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([55.0,147.0],
                   [147.0,147.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([55.0,147.0],
                   [55.0,55.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([162.0,162.0],
                   [147.0,162.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([147.0,147.0],
                   [147.0,162.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([147.0,162.0],
                   [162.0,162.0],
                   color='black',
                   ls=':',
                   lw=0.5)
    PylabAxis.plot([147.0,162.0],
                   [147.0,147.0],
                   color='black',
                   ls=':',
                   lw=0.5)
#    PylabImAxis.set_ylim(PylabImAxis.get_ylim()[::-1])
#    PylabAxis.set_ylim(PylabAxis.get_ylim()[::-1])
    PylabImAxis.set_xlim(Extent[0:2])
    PylabAxis.set_xlim(Extent[0:2])
    PylabImAxis.set_ylim(Extent[2:4])
    PylabAxis.set_ylim(Extent[2:4])

    PylabColorBar = pylab.colorbar(PylabImage)
    PylabColorBar.set_label(ZLabel)

    pylab.savefig('HeatMap.pdf')
    pylab.show()

    return

def main(ExecutableName):

    ArgParser,\
    Arguments   = ParseArguments()

    Ext = ''
    Log = Logger.Logger(ExecutableName,
                        Ext)
    LogString  = '## START TIMESTAMP\n'
    LogString += str(Log.GetStartDate())+'\n'
    LogString += '## END TIMESTAMP'
    print LogString
    Log.Write(LogString+'\n')
    LogString = Log.GetStartLogString()
    print LogString
    Log.Write(LogString+'\n')
    LogArguments(Log,
                 ArgParser,
                 Arguments)

    LogString  = '**** Parsing correlation data file \"'+Arguments.RatioFileName+'\" ...\n'
    LogString += '**** Writing sorted data to \"'+re.sub('.txt','.sorted.txt',Arguments.RatioFileName)+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    fr   = open(Arguments.RatioFileName,'r')

    IndexA     = []
    IndexB     = []
    NameA      = []
    NameB      = []
    Slopes     = []
    Intercepts = []
    R2s        = []
    PVals      = []
    StdErrs    = []
    fw         = open(re.sub('.txt','.sorted.txt',Arguments.RatioFileName),'w')
    for Line in fr:
        if(Line[0]!='#'):
            LSplit = Line.strip().split()
            IndexA.append(int(LSplit[0]))
            IndexB.append(int(LSplit[1]))
            NameA.append(LSplit[2])
            NameB.append(LSplit[3])
            Slopes.append(float(LSplit[4]))
            Intercepts.append(float(LSplit[5]))
            R2s.append(float(LSplit[6]))
            PVals.append(float(LSplit[7]))
            StdErrs.append(float(LSplit[8]))
        else:
            fw.write(Line)
    fr.close()

    SortIndices = list(scipy.argsort(scipy.array(R2s)))
    SortIndices.reverse() # highest R2 values first
    for i in SortIndices:
        FormatString = '{0:>8s}'
        fw.write(FormatString.format(str(IndexA[i])))
        FormatString = '{0:>7s}'
        fw.write(FormatString.format(str(IndexB[i])))
        FormatString = '{0:>16s}'
        fw.write(FormatString.format(NameA[i]))
        fw.write(FormatString.format(NameB[i]))
        FormatString = '{0:>+14.5e}'
        fw.write(FormatString.format(Slopes[i]))
        fw.write(FormatString.format(Intercepts[i]))
        fw.write(FormatString.format(R2s[i]))
        fw.write(FormatString.format(PVals[i]))
        fw.write(FormatString.format(StdErrs[i]))
        fw.write('\n')
    fw.close()

    X = scipy.array(IndexA)
    Y = scipy.array(IndexB)
    Z = scipy.array(R2s)*scipy.sign(scipy.array(Slopes))
#    Z = scipy.array(R2s)*scipy.sign(scipy.array(Slopes))
    PlotHeatMap(X=X,
                Y=Y,
                Z=Z,
                XLabel=r'${\rm Index~Metabolite~A}$',
                YLabel=r'${\rm Index~Metabolite~B}$',
                ZLabel=r'${\rm sign}(a)R^{2}$')

    LogString = '**** Done :-)'
    print LogString
    Log.Write(LogString+'\n')
    LogString = Log.GetEndLogString()
    print LogString
    Log.Write(LogString+'\n')
    Log.Close()

    return

if(__name__=='__main__'):
    ExecutableName = os.path.abspath(__file__).split('/')[-1]
    main(ExecutableName)
