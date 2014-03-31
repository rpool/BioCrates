#! /usr/bin/env python
import os
import re
import colorsys
import scipy
import sys
import pylab
import argparse

import ArgumentParser
import Logger
import DataContainer

def PylabUpdateParams():
    figwidth_pt   = 1422.0 # pt (from revtex \showthe\columnwidth)
#        figwidth_pt   = 246.0 # pt (from revtex \showthe\columnwidth)
    inches_per_pt = 1.0/72.27
    figwidth      = figwidth_pt*inches_per_pt
    golden_mean   = (scipy.sqrt(5.0)-1.0)/2.0 # Aesthetic ratio
    figheight     = figwidth*golden_mean
    fig_size      = [figwidth,figheight]
    params        = {'backend': 'pdf',
                     'patch.antialiased': True,
                     'axes.labelsize': 28,
                     'axes.linewidth': 0.5,
                     'grid.color': '0.75',
                     'grid.linewidth': 0.25,
                     'grid.linestyle': ':',
                     'axes.axisbelow': False,
                     'text.fontsize': 24,
                     'legend.fontsize': 24,
                     'xtick.labelsize': 24,
                     'ytick.labelsize': 24,
                     'text.usetex': True,
                     'figure.figsize': fig_size}
    left   = 0.06
    bottom = 0.125
    width  = 0.95-left
    height = 0.95-bottom

    return params,\
           [left, bottom, width, height]

def PlotManhattan(xname=str,
                  yname=str,
                  DCs=DataContainer.DataContainers,
                  Arguments=argparse.Namespace,
                  PhenotypeName=str,
                  SNPInfoDict=dict,
                  XPosArray=scipy.array,
                  ChrArray=scipy.array,
                  XTicks=list,
                  XTickLabels=list,
                  XXMin=int,
                  XXMax=int,
                  Log=Logger):

    LogString = '**** Generating Manhattan plot ...'
    print LogString
    Log.Write(LogString+'\n')

    PylabParameters,\
    Rectangle         = PylabUpdateParams()
    pylab.rcParams.update(PylabParameters)
    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis   = PylabFigure.add_axes(Rectangle)
#    XMax        = 0
#    XTicks      = []
#    XTickLabels = []
#    for i in range(Arguments.NChr):
#        XName = ''
#        YName = ''
#    for Key in DCs.DataContainers.iterkeys():
#        if(re.search(xname,Key)):
#            XName = Key
#        if(re.search(yname,Key)):
#            YName = Key
    SNPIDArray = DCs.DataContainers['MarkerName'].GetDataArray()
    PValArray  = scipy.array(DCs.DataContainers['P-value'].GetDataArray()).astype(float)
    Y          = -scipy.log10(PValArray)

    IndexArray = []
    for Entry in SNPIDArray:
        IndexArray.append(SNPInfoDict[Entry]['index'])
    IndexArray = scipy.array(IndexArray)
    X          = XPosArray[IndexArray]
    ChromArray = ChrArray[IndexArray]
    Colors = ['black',
              'grey',
              'black',
              'grey',
              'black',
              'grey',
              'black',
              'grey',
              'black',
              'grey',
              'black',
              'grey',
              'black',
              'grey',
              'black',
              'grey',
              'black',
              'grey',
              'black',
              'grey',
              'black',
              'grey']
    ColorArray = []
    for Entry in ChromArray:
        ColorArray.append(Colors[int(Entry)-1])
    ColorArray = scipy.array(ColorArray)

    YInsign = Y < -scipy.log10(1.0e-6/13.0)
    YSign   = Y > -scipy.log10(5.0e-8/13.0)
    YSugg   = Y >= -scipy.log10(1.0e-6/13.0)
    YSugg  *= Y <= -scipy.log10(5.0e-8/13.0)

    YY = scipy.compress(YInsign,Y)
    if(len(YY)>0):
        PylabAxis.scatter(x=scipy.compress(YInsign,X),
                          y=YY,
                          color=scipy.compress(YInsign,ColorArray).tolist(),
                          s=0.5)
    YY = scipy.compress(YSugg,Y)
    if(len(YY)>0):
        PylabAxis.scatter(x=scipy.compress(YSugg,X),
                          y=YY,
                          color=scipy.compress(YSugg,ColorArray).tolist(),
                          s=5.0)
    YY = scipy.compress(YSign,Y)
    if(len(YY)>0):
        PylabAxis.scatter(x=scipy.compress(YSign,X),
                          y=YY,
                          color=scipy.compress(YSign,ColorArray).tolist(),
                          s=10.0)

    XSign = scipy.array(PylabAxis.get_xlim())
    YSign = -scipy.log10(scipy.array([5.0e-8/13.0,5.0e-8/13.0]))
    PylabAxis.plot(XSign,
                   YSign,
                   linestyle='--',
                   color='grey',
                   label=r'${\rm '+PhenotypeName+'}$',
                   linewidth=1.25)
    XSugg = scipy.array(PylabAxis.get_xlim())
    YSugg = -scipy.log10(scipy.array([1.0e-6/13.0,1.0e-6/13.0]))
    PylabAxis.plot(XSugg,
                   YSugg,
                   linestyle=':',
                   color='grey',
                   linewidth=1.25)
#    PylabAxis.set_ylim([0.0,PylabAxis.get_ylim()[1]])
#    PylabAxis.set_xlim([0,XMax])
    Handles,Labels = PylabAxis.get_legend_handles_labels()
    PylabAxis.legend(Handles,
                     Labels,
                     fancybox=True,
                     shadow=True,
                     loc='best')
    PylabAxis.set_xlabel(r'$\rm position$')
    PylabAxis.spines['right'].set_visible(False)
    PylabAxis.spines['top'].set_visible(False)
    PylabAxis.xaxis.set_ticks_position('bottom')
    PylabAxis.yaxis.set_ticks_position('left')
    PylabAxis.set_ylim([0.0,PylabAxis.get_ylim()[1]])
    XXRange  = float(XXMax)-float(XXMin)
    XXOffset = XXRange*0.005
    PylabAxis.set_xlim([float(XXMin)-XXOffset,float(XXMax)+XXOffset])
    PylabAxis.xaxis.set_ticks(XTicks)
    PylabAxis.xaxis.set_ticklabels(XTickLabels)
    for Label in PylabAxis.xaxis.get_ticklabels():
        Label.set_rotation(90)
#    PylabAxis.xaxis.set_ticks(XTicks)
#    PylabAxis.xaxis.set_ticklabels(XTickLabels)
    for Label in PylabAxis.xaxis.get_ticklabels():
        Label.set_rotation(90)
    PylabAxis.set_ylabel(r'$-{\rm log}_{10}(p-{\rm value})$')

    PylabFigure.savefig('Manhattan_'+PhenotypeName+'.png')
    PylabAxis.clear()
    pylab.close(PylabFigure)
    del PylabFigure
    del PylabAxis

    return

def main(ExecutableName):

    ArgParser,\
    Arguments   = ArgumentParser.ParseArguments()

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
    ArgumentParser.LogArguments(Log,
                                ArgParser,
                                Arguments)

    TmpPath = os.path.join(os.getcwd(),'Tmp')
    if(not os.path.isdir(TmpPath)):
        os.mkdir(TmpPath)
    SNPInfoDict  = {} # will contain a dictionary of two elements: {'chr':int(chr),'pos':int(pos)}
    SNPChrArray  = scipy.array([])
    SNPPosArray  = scipy.array([])
    SNPXPosArray = scipy.array([])
    XLeft        = []
    XRight       = []
    XTicks       = []
    XTickLabels  = []
    XXMax        = None
    XXMin        = None
    XPath               = os.path.join(Arguments.GWAOutputPath,'SNPInfo.txt')
    SNPInfoDecompressed = os.path.join(TmpPath,'SNPInfo.txt')
    os.system('ln -sf '+XPath+' '+SNPInfoDecompressed)
#        os.system('pigz -d -c -k '+XPath+' > '+SNPInfoDecompressed)
    FH          = open(SNPInfoDecompressed,'r')
    HeaderList  = FH.readline().strip().split()
    FH.close()
    SNPIDArray  = None
    ChrArray    = None
    PosArray    = None
    SNPIDColumn = None
    ChrColumn   = None
    PosColumn   = None
    for Entry in HeaderList:
        if(Entry=='SNPID'):
            SNPIDColumn = HeaderList.index(Entry)
        if(Entry=='chr'):
            ChrColumn = HeaderList.index(Entry)
        if(Entry=='position'):
            PosColumn = HeaderList.index(Entry)
    Arrays = scipy.loadtxt(fname=SNPInfoDecompressed,
                           dtype=str,
                           skiprows=1,
                           usecols=[SNPIDColumn,
                                    ChrColumn,
                                    PosColumn],
                           unpack=True)
    os.remove(SNPInfoDecompressed)
    SNPIDArray = Arrays[0]
    ChrArray   = Arrays[1].astype(int)
    PosArray   = Arrays[2].astype(int)
    ArgSort    = scipy.argsort(ChrArray)
    SNPIDArray = SNPIDArray[ArgSort]
    ChrArray   = ChrArray[ArgSort]
    PosArray   = PosArray[ArgSort]

    TmpSNPIDArray = scipy.array([])
    TmpChrArray   = scipy.array([])
    TmpPosArray   = scipy.array([])
    TmpXPosArray  = scipy.array([])
    XMax          = 0
    for c in range(Arguments.NChr):
        Chr              = c+1
        FilterArray      = (ChrArray==Chr)
        TmpTmpSNPIDArray = scipy.compress(FilterArray,SNPIDArray)
        TmpTmpChrArray   = scipy.compress(FilterArray,ChrArray)
        TmpTmpPosArray   = scipy.compress(FilterArray,PosArray)
        ArgSort          = scipy.argsort(TmpTmpPosArray)
        TmpTmpSNPIDArray = TmpTmpSNPIDArray[ArgSort]
        TmpTmpChrArray   = TmpTmpChrArray[ArgSort]
        TmpTmpPosArray   = TmpTmpPosArray[ArgSort]
        TmpTmpXPosArray  = TmpTmpPosArray+XMax
        XLeft.append(float(TmpTmpXPosArray.min()))
        XRight.append(float(TmpTmpXPosArray.max()))
        XTicks.append(0.5*(XLeft[-1]+XRight[-1]))
        XTickLabels.append(r'${\rm CHR'+str(Chr)+r'}$')

        TmpSNPIDArray = scipy.append(TmpSNPIDArray,TmpTmpSNPIDArray)
        TmpChrArray   = scipy.append(TmpChrArray,TmpTmpChrArray)
        TmpPosArray   = scipy.append(TmpPosArray,TmpTmpPosArray)
        TmpXPosArray  = scipy.append(TmpXPosArray,TmpTmpXPosArray)
        XMax          = TmpXPosArray.max()
        del TmpTmpSNPIDArray
        del TmpTmpChrArray
        del TmpTmpPosArray
        del TmpTmpXPosArray

    SNPIDArray = scipy.array(TmpSNPIDArray)
    ChrArray   = scipy.array(TmpChrArray)
    PosArray   = scipy.array(TmpPosArray)
    XPosArray  = scipy.array(TmpXPosArray)
    del TmpSNPIDArray
    del TmpChrArray
    del TmpPosArray
    del TmpXPosArray

    for i in range(len(SNPIDArray)):
        Entry                       = SNPIDArray[i]
        SNPInfoDict[Entry]          = {}
        SNPInfoDict[Entry]['index'] = i
    SNPChrArray  = scipy.array([])
    SNPPosArray  = scipy.array([])
    SNPXPosArray = scipy.array([])
    SNPChrArray  = scipy.append(SNPChrArray,ChrArray).astype(int)
    SNPPosArray  = scipy.append(SNPPosArray,PosArray).astype(int)
    SNPXPosArray = scipy.append(SNPXPosArray,XPosArray).astype(int)
    XXMin = XPosArray.min()
    XXMax = XPosArray.max()
    LogString = '**** Parsed \"SNPInfo.txt\" ...'
    print LogString
    Log.Write(LogString+'\n')

    if(Arguments.GWAOutputFile!=None):
        DCs = DataContainer.DataContainers()
        DCs.ParseGWAOutput(Arguments.SnpTestOutputFile,
                           Log)
        DCs.PlotManhattan(xname='pos',
                          yname='pvalue',
                          Log=Log)
    else:
        GWAOutputFiles = []
        for File in os.listdir(Arguments.GWAOutputPath):
            if((re.search(Arguments.GWAOutputExtStr,File)) and
               (not re.search('.log',File))):
                GWAOutputFiles.append(File)
        for p in range(len(GWAOutputFiles)):
            File  = GWAOutputFiles[p]
            P     = GWAOutputFiles[p].split('_')[1]

#            DCsList = DataContainer.ListDataContainers()
#            DCsList.SetPhenotypeName(P)

            HSV_tuples = None
            RGB_tuples = None
            if(Arguments.boGreyScale):
                boGrey     = True
                GreyHSV    = (0.0,0.0,0.4)
                BlackHSV   = (0.0,0.0,0.0)
                HSV_tuples = []
                for i in range(Arguments.NChr):
                    if(boGrey):
                        HSV_tuples.append(GreyHSV)
                        boGrey = False
                    else:
                        HSV_tuples.append(BlackHSV)
                        boGrey = True
                RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
            else:
                HSV_tuples = [(x*1.0/Arguments.Nchr, 0.75, 0.75) for x in range(Arguments.NChr)]
                RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)

            DCs = DataContainer.DataContainers()
            DCs.ParseGWAOutput(os.path.join(Arguments.GWAOutputPath,File),
                               Log)
            DCs.Color = RGB_tuples[i]
#            DCsList.List.append(DCs)
            PlotManhattan('pos',
                          'P-value',
                          DCs,
                          Arguments,
                          P,
                          SNPInfoDict,
                          XPosArray,
                          ChrArray,
                          XTicks,
                          XTickLabels,
                          XXMin,
                          XXMax,

                          Log)
#            DCsList.PlotManhattan(xname='pos',
#                                  yname='P-value',
#                                  Log=Log)

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