#! /usr/bin/env python

# python modules
import os
import fnmatch
import sys
import scipy
import scipy.linalg
import scipy.interpolate
import re
import matplotlib.mlab
import pylab

# Homebrew modules
import Logger
import ArgumentParser
import File
import Plotting

def PylabGetParams():
#   Set the plotting parameters: TeX mode, two columns per page in revtex mode/
    FigwidthPt  = 246.0 # pt (from revtex \showthe\columnwidth)
    InchesPerPt = 1.0/72.27
    FigWidth    = FigwidthPt*InchesPerPt
    GoldenMean  = (scipy.sqrt(5.0)-1.0)/2.0 # Aesthetic ratio
    FigHeight   = FigWidth*GoldenMean
    FigSize     = [FigWidth,FigHeight]
    Params      = {'backend': 'pdf',
                   'patch.antialiased': True,
                   'axes.labelsize': 8,
                   'axes.linewidth': 0.5,
                   'axes.fontsize' : 4,
                   'grid.color': '0.75',
                   'grid.linewidth': 0.125,
                   'grid.linestyle': ':',
                   'axes.axisbelow': False,
                   'text.fontsize': 8,
                   'legend.fontsize': 4,
                   'xtick.labelsize': 8,
                   'ytick.labelsize': 8,
                   'text.usetex': True,
                   'figure.figsize': FigSize}
    Left   = 0.15
    Bottom = 0.16
    Width  = 0.84 - Left
    Height = 0.95 - Bottom

    return Params,\
           [Left, Bottom, Width, Height]

def ScreeTest(Data=scipy.array):
    U, S, V     = scipy.linalg.svd(Data)
    EigenValues = scipy.power(S,2.0) / scipy.sum(S)

    PylabParameters,\
    Rectangle         = PylabGetParams()
    Size              = 1.5
    LineWidth         = 0.5
    pylab.rcParams.update(PylabParameters)

    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis = PylabFigure.add_axes(Rectangle)

    SingVals = scipy.arange(len(EigenValues)) + 1
    PylabAxis.plot(SingVals,
                   EigenValues,
                   'ro-',
                   linewidth=0.75,
                   markersize=2.0)

    for i in range(len(SingVals)):
        if(EigenValues[i]<1.0):
            KaiserCriterion = SingVals[i]
            break

    dX   = scipy.diff(SingVals)
    dY   = scipy.diff(EigenValues)
    dYdX = dY/dX
#    PylabAxis.plot(SingVals[1:],
#                   dYdX,
#                   'g-',
#                   linewidth=0.75,
#                   markersize=2.0)
    dX2    = scipy.diff(dX)
    dY2    = scipy.diff(dY)
    dY2dX2 = dY2/dX[1:]
    PylabAxis.plot(SingVals[2:],
                   dY2dX2,
                   'b-',
                   linewidth=0.75,
                   markersize=2.0)
    ElbowPoint = None
    for i in range(2,len(SingVals)):
        if(dY2dX2[i]<0.0):
            ElbowPoint = SingVals[i]
            break


    PylabAxis.plot(scipy.array([0.0,SingVals.max()]),
                   scipy.array([1.0,1.0]),
                   '-',
                   color='grey',
                   linewidth=0.5,
                   markersize=2.0)
    PylabAxis.set_xlabel(r'$\rm Principal~Component$')
    PylabAxis.set_ylabel(r'$\rm Eigenvalue$')
    Legend = pylab.legend([r'$\rm Eigenvalues~from~SVD$',
                           r'$\rm Second~derivative$'],
                           loc='best',
                           fancybox=True,
                           shadow=True)
    PylabAxis.grid(True)
    pylab.savefig('ScreePlot.png',
                  dpi=600)
#    pylab.show()

    return ElbowPoint,\
           KaiserCriterion

def main(ExecutableName=str):

    #===========================================================================
    # START Initialization
    #===========================================================================

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
    #===========================================================================
    # END Initialization
    #===========================================================================

    #===========================================================================
    # START Do the work!
    #===========================================================================

    # Parse Arguments.ProtocolFile
    LogString = '++ Parsing \"'+Arguments.ProtocolFile+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    ProtocolFile = File.File(Name=Arguments.ProtocolFile,
                             boHeader=False)
    XmlProtocol  = ProtocolFile.ParseAsXml()
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    ExcludedMtbs = []
    for Tag in XmlProtocol.getroot().find('ExcludedMtbs'):
        ExcludedMtbs.append(Tag.text)

    PhenotypePath       = XmlProtocol.getroot().find('PhenotypeInputPath').text.strip()
    PhenotypeExtension  = XmlProtocol.getroot().find('PhenotypeFileExt').text.strip()
    PhenotypeInputFiles = fnmatch.filter(os.listdir(PhenotypePath),'*.'+PhenotypeExtension)
    PhenotypeInputFiles = fnmatch.filter(PhenotypeInputFiles,'*Phenotype*'+PhenotypeExtension)
    fr = open(PhenotypeInputFiles[0],'r')
    HeaderList     = fr.readline().strip().split(',')
    MetaboliteList = []
    fr.close()
    for i in range(len(HeaderList)):
        Entry         = re.sub('\"','',HeaderList[i])
        Entry         = re.sub('\.\.','.',Entry)
        HeaderList[i] = Entry
        if(i>1):
            MetaboliteList.append(Entry)
    UseCols = range(1,len(HeaderList))
    Arrays = scipy.loadtxt(fname=PhenotypeInputFiles[0],
                           dtype=str,
                           delimiter=',',
                           skiprows=1,
                           usecols=UseCols,
                           unpack=True)
    SampleIdArray = Arrays[0].astype(int)

    PhenotypeArrayDict = {}
    RawData            = []
    for i in range(len(MetaboliteList)):
        Mtb = MetaboliteList[i]
        PhenotypeArrayDict[Mtb] = Arrays[i+1].astype(float)
        RawData.append(PhenotypeArrayDict[Mtb])

    MeanCenteredData           = []
    MeanCenteredAutoScaledData = []
    for i in range(len(RawData)):
        Array = RawData[i]
        Mean  = scipy.mean(Array)
        Std   = scipy.std(Array,ddof=1)
        MeanCenteredData.append(Array-Mean)
        MeanCenteredAutoScaledData.append(MeanCenteredData[-1]/Std)

    LogString = '++ Performing PCA analysis on raw data array ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString = '  ++ Writing fraction explained variance results to \"RawPCA.out.csv\" ...'
    print LogString
    Log.Write(LogString+'\n')
    fw = open('RawPCA.out.csv','w')
    fw.write('PC,Frac,CummFrac\n')
    MyData = scipy.array(RawData).T
    MyPCA  = matplotlib.mlab.PCA(MyData)
    Sum    = 0.0
    for i in range(len(MyPCA.fracs)):
        Sum += MyPCA.fracs[i]
        fw.write(str(i+1)+','+str(MyPCA.fracs[i])+','+str(Sum)+'\n')
    fw.close()
    LogString = '  -- Done ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    LogString = '++ Performing PCA analysis on mean centered data array ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString = '  ++ Writing fraction explained variance results to \"MeanCenteredPCA.out.csv\" ...'
    print LogString
    Log.Write(LogString+'\n')
    fw = open('MeanCenteredPCA.out.csv','w')
    fw.write('PC,Frac,CummFrac\n')
    MyData = scipy.array(MeanCenteredData).T
    MyPCA  = matplotlib.mlab.PCA(MyData)
    Sum    = 0.0
    for i in range(len(MyPCA.fracs)):
        Sum += MyPCA.fracs[i]
        fw.write(str(i+1)+','+str(MyPCA.fracs[i])+','+str(Sum)+'\n')
    fw.close()
    LogString = '  -- Done ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    LogString = '++ Performing PCA analysis on mean centered auto scaled data array ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString = '  ++ Writing fraction explained variance results to \"MeanCenteredAutoScaledPCA.out.csv\" ...'
    print LogString
    Log.Write(LogString+'\n')
    fw = open('MeanCenteredAutoScaledPCA.out.csv','w')
    fw.write('PC,Frac,CummFrac\n')
    MyData = scipy.array(MeanCenteredAutoScaledData).T
    MyPCA  = matplotlib.mlab.PCA(MyData)
    Sum    = 0.0
    for i in range(len(MyPCA.fracs)):
        Sum += MyPCA.fracs[i]
        fw.write(str(i+1)+','+str(MyPCA.fracs[i])+','+str(Sum)+'\n')
    fw.close()
    LogString = '  -- Done ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    LogString = '++ Performing scree test and determining Kaiser criterion ...'
    print LogString
    Log.Write(LogString+'\n')
    ElbowPoint,\
    KaiserCriterion = ScreeTest(MeanCenteredAutoScaledData)
    LogString = '  ** Saved plot to \"ScreePlot.png\" ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString  = '  ** Elbow point in scree plot = '+str(ElbowPoint)+' ...\n'
    LogString += '  ** This is defined as the point after which the numerical second derivate of the eigenvalues(#PC) goes below 0.0 ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString  = '  ** Kaiser critrion = '+str(KaiserCriterion)+' ...\n'
    LogString += '  ** This is defined as the point after which eigenvalues(#PC) goes below 1.0 ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    #===========================================================================
    # END Do the work!
    #===========================================================================

    #===========================================================================
    # START Finalization
    #===========================================================================
    LogString = '\n**** Done :-)'
    print LogString
    Log.Write(LogString+'\n')
    LogString = Log.GetEndLogString()
    print LogString
    Log.Write(LogString+'\n')
    Log.Close()
    #===========================================================================
    # START Finalization
    #===========================================================================

    return

if(__name__=='__main__'):
    ExecutableName = os.path.abspath(__file__).split('/')[-1]
    main(ExecutableName)