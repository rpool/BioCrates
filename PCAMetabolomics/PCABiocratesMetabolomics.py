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

def PlotScoresAndLoadings(Data=scipy.array,
                          NPCs=int,
                          PCA=matplotlib.mlab.PCA,
                          LoadingPlotColors=list,
                          LoadingLegendList=list,
                          MtbName2ClassesDict=dict,
                          MtbList=list):

    U, S, V     = scipy.linalg.svd(Data,full_matrices=False)
    EigenValues = scipy.power(S,2.0)/float(len(U))

    IndexesList = []
    for i in range(len(LoadingLegendList)):
        IndexesList.append([])
    for i in range(len(MtbList)):
        LoadingLegendListIndex = LoadingLegendList.index(MtbName2ClassesDict[MtbList[i]])
        IndexesList[LoadingLegendListIndex].append(i)

    for i in range(NPCs-1):
        for j in range(i+1,NPCs):
            pylab.close()

            PylabParameters,\
            Rectangle         = PylabGetParams()
            Size              = 1.5
            LineWidth         = 0.5
            pylab.rcParams.update(PylabParameters)

            Fig = pylab.figure(figsize=(14.0*0.394,7.0*0.394))
            Fig.clf()
            Ax1 = Fig.add_subplot(121)
            Ax2 = Fig.add_subplot(122)

            Ax1.scatter(PCA.Y[:,i],
                        PCA.Y[:,j],
                        s=Size,
                        alpha=0.75,
                        facecolor='None',
                        linewidth=0.75,
                        color='grey')

            XLine  = scipy.array([PCA.Y[:,i].min(),PCA.Y[:,i].max()])
            XLine /= PCA.Y[:,i].max()-PCA.Y[:,i].min()
            XLine *= EigenValues[i]
            YLine = scipy.array([PCA.Y[:,j].min(),PCA.Y[:,j].max()])
            YLine /= PCA.Y[:,j].max()-PCA.Y[:,j].min()
            YLine *= EigenValues[j]
            Ax1.plot(XLine,
                     scipy.array([0.0,0.0]),
                     'k-',
                     lw=LineWidth)
            Ax1.plot(scipy.array([0.0,0.0]),
                     YLine,
                     'k-',
                     lw=LineWidth)
            Ax1.set_xlabel(r'$\rm Scores~PC'+str(i+1)+r'~('+str(round(PCA.fracs[i]*100.0,2))+r'\%)$',size=4)
            Ax1.set_ylabel(r'$\rm Scores~PC'+str(j+1)+r'~('+str(round(PCA.fracs[j]*100.0,2))+r'\%)$',size=4)
            for tick in Ax1.xaxis.get_major_ticks():
                tick.label.set_fontsize(4)
            for tick in Ax1.yaxis.get_major_ticks():
                tick.label.set_fontsize(4)

            Ax1.grid(True)
            Ax1.set_axisbelow(True)

            for I in range(len(LoadingLegendList)):
                Ax2.scatter(PCA.Wt[i,scipy.array(IndexesList[I])],
                            PCA.Wt[j,scipy.array(IndexesList[I])],
                            color=LoadingPlotColors[I],
                            s=Size,
                            facecolor='None',
                            linewidth=0.75,
                            alpha=0.75)
        #    for i in range(len(LoadingLegendList)-1)
        #        Ax2.scatter(PCA.Wt[:,0],
        #                    PCA.Wt[:,1],
        #                    color=LoadingPlotColors,
        #                    s=Size,
        #                    alpha=0.5) # dummy plots for legend

            LegendList = []
            for Entry in LoadingLegendList:
                LegendList.append(r'$\rm '+r'~'.join(Entry.split())+r'$')

            Ax2.legend(LegendList,
                       loc='best',
                       fancybox=True,
                       shadow=True,
                       scatterpoints=1)

            Ax2.plot(Ax2.get_xlim(),
                     scipy.array([0.0,0.0]),
                     '-',
                     color='grey',
                     lw=LineWidth,
                     zorder=0)
            Ax2.plot(scipy.array([0.0,0.0]),
                     Ax2.get_ylim(),
                     '-',
                     color='grey',
                     lw=LineWidth,
                     zorder=0)

            Ax2.set_xlabel(r'$\rm Loadings~PC'+str(i+1)+r'~('+str(round(PCA.fracs[i]*100.0,2))+r'\%)$',size=4)
            Ax2.set_ylabel(r'$\rm Loadings~PC'+str(j+1)+r'~('+str(round(PCA.fracs[j]*100.0,2))+r'\%)$',size=4)
            for tick in Ax2.xaxis.get_major_ticks():
                tick.label.set_fontsize(4)
            for tick in Ax2.yaxis.get_major_ticks():
                tick.label.set_fontsize(4)

            Ax2.grid(True)
            Ax2.set_axisbelow(True)

            PlotName = 'ScoresAndLoadingsPC'+str(i+1)+'_PC'+str(j+1)+'.png'
            Fig.subplots_adjust(wspace=0.4)
            Fig.savefig(PlotName,dpi=1000)

    return

def HornsParallalAnalysis(DataDict=dict,
                          Log=Logger):
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    import rpy2.rinterface as rinterface

    Paran = None

    LogString  = '  ** Writing rpy2 parallel analysis output to \"'+Log.GetFileName()+'\" ...\n'
    LogString += '  ## START rpy2 ##'
    print LogString
    Log.Write(LogString+'\n')

    StdOutSav     = sys.stdout
    sys.stdout    = Log.GetFileHandle()
    DataFrameDict = {}
    for Key,Value in DataDict.iteritems():
        hlen = len(Value)
        DataFrameDict[Key] = rinterface.baseenv['as.real'](rinterface.StrSexpVector(Value[0:hlen]))
        RDataFrame  = robjects.DataFrame(DataFrameDict)
    Paran       = importr('paran')
#    ParanOutput = Paran.paran(RDataFrame,iterations=180) # for speed
    ParanOutput = Paran.paran(RDataFrame) # default number of iterations

    sys.stdout  = StdOutSav
    LogString   = '  ## END rpy2 ##'
    print LogString
    Log.Write(LogString+'\n')

    LogString  = '  ** Plotting parallel analysis output to \"ParanPlot.png\" ...'
    print LogString
    Log.Write(LogString+'\n')
    PylabParameters,\
    Rectangle         = PylabGetParams()
    Size              = 1.5
    LineWidth         = 0.5
    pylab.rcParams.update(PylabParameters)

    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis = PylabFigure.add_axes(Rectangle)

    NRetained = int(scipy.array(ParanOutput[ParanOutput.names.index('Retained')])[0])

    Y = scipy.array(ParanOutput[ParanOutput.names.index('RndEv')])
    X = scipy.arange(len(Y))+1
    PylabAxis.plot(X,
                   Y,
                   'bs-',
                   markeredgecolor='blue',
                   linewidth=0.5,
                   alpha=0.5,
                   markersize=2.0)

    Y = scipy.array(ParanOutput[ParanOutput.names.index('Ev')])
    X = scipy.arange(len(Y))+1
    PylabAxis.plot(X,
                   Y,
                   'ro-',
                   markeredgecolor='red',
                   linewidth=0.5,
                   alpha=0.5,
                   markersize=2.0)

    Y = scipy.array(ParanOutput[ParanOutput.names.index('AdjEv')])
    X = scipy.arange(len(Y))+1
    PylabAxis.plot(X[:NRetained],
                   Y[:NRetained],
                   'ko-',
                   linewidth=0.5,
                   markersize=2.0)
    PylabAxis.plot(X[NRetained:-1],
                   Y[NRetained:-1],
                   'ko-',
                   markerfacecolor='None',
                   linewidth=0.5,
                   markersize=2.0)

    PylabAxis.plot(scipy.array([0.0,X.max()]),
                   scipy.array([1.0,1.0]),
                   '-',
                   color='grey',
                   linewidth=0.5)

    PylabAxis.set_xlabel(r'$\rm Principal~Component$')
    PylabAxis.set_ylabel(r'$\rm Eigenvalue$')
    Legend = pylab.legend([r'$\rm Random~eigenvalues$',
                           r'$\rm Unadjusted~eigenvalues$',
                           r'$\rm Adjusted~eigenvalues~(retained~PCs)$',
                           r'$\rm Adjusted~eigenvalues~(unretained~PCs)$'],
                           loc='best',
                           fancybox=True,
                           shadow=True)
    PylabAxis.grid(True)
    pylab.savefig('ParanPlot.png',
                  dpi=600)

    return NRetained

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
    EigenValues = scipy.power(S,2.0)/float(len(U))

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
    PhenotypeBName      = XmlProtocol.getroot().find('PhenotypeFileBase').text.strip()
    PhenotypeExtension  = XmlProtocol.getroot().find('PhenotypeFileExt').text.strip()
    PhenotypeInputFile  = PhenotypeBName+'.'+PhenotypeExtension
    fr = open(PhenotypeInputFile,'r')
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
    LogString = '++ Parsing \"'+PhenotypeInputFile+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    Arrays = scipy.loadtxt(fname=PhenotypeInputFile,
                           dtype=str,
                           delimiter=',',
                           skiprows=1,
                           usecols=UseCols,
                           unpack=True)
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')
    SampleIdArray = Arrays[0].astype(int)

    MtbName2ClassesDict = {}
    MtbClasses          = []
    if(eval(XmlProtocol.getroot().find('MtbClassesFile').find('boUse').text)):
        FName = XmlProtocol.getroot().find('MtbClassesFile').find('Path').text
        FName = os.path.join(FName,XmlProtocol.getroot().find('MtbClassesFile').find('Name').text)
        fr = open(FName,'r')
        for Line in fr:
            if(Line[0]=='#'):
                continue
            else:
                MtbName = Line.strip().split()[1]
                MtbClass = Line.strip().split()[2]
                MtbName2ClassesDict[MtbName] = MtbClass
                MtbClasses.append(MtbClass)
        fr.close()
        MtbClasses = list(set(MtbClasses))
    MtbClassColors = ['black',
                      'red',
                      'blue',
                      'green',
                      'purple']

    PhenotypeArrayDict = {}
    RawData            = []
    for i in range(len(MetaboliteList)):
        Mtb = MetaboliteList[i]
        PhenotypeArrayDict[Mtb] = Arrays[i+1]
        RawData.append(PhenotypeArrayDict[Mtb])

    LogString = '++ Accounting for excluded metabolites in current data set ...'
    print LogString
    Log.Write(LogString+'\n')
    DelList = []
    for Key in PhenotypeArrayDict.iterkeys():
        if(Key in ExcludedMtbs):
            LogString = '  ** Metabolite \"'+Key+'\" will be excluded from the data set ...'
            print LogString
            Log.Write(LogString+'\n')
            DelList.append(Key)
    if(len(DelList)==0):
        LogString = '  ** No metabolites are exluded ...'
        print LogString
        Log.Write(LogString+'\n')
    else:
        LogString = '  ** '+str(len(DelList))+' metabolites are exluded ...'
        print LogString
        Log.Write(LogString+'\n')

    for Entry in DelList:
        del PhenotypeArrayDict[Entry]
    del DelList
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    RawData            = []
    MetaboliteList     = PhenotypeArrayDict.keys()
    for i in range(len(MetaboliteList)):
        Mtb = MetaboliteList[i]
        RawData.append(PhenotypeArrayDict[Mtb].astype(float))

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
    KaiserCriterion = ScreeTest(MyData)
    LogString = '  ** Saved plot to \"ScreePlot.png\" ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString  = '  ** Elbow point in scree plot = '+str(ElbowPoint)+' ...\n'
    LogString += '  ** This is defined as the point after which the numerical second derivate of the eigenvalues(#PC) goes below 0.0 ...\n'
    LogString += '  ** The percentage explained variance for this number of PCs is '+str(round(scipy.sum(MyPCA.fracs[:ElbowPoint])*100.0,2))+'% (based upon mean centered auto scaled input data) ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString  = '  ** Kaiser critrion = '+str(KaiserCriterion)+' ...\n'
    LogString += '  ** This is defined as the point after which eigenvalues(#PC) goes below 1.0 ...\n'
    LogString += '  ** The percentage explained variance for this number of PCs is '+str(round(scipy.sum(MyPCA.fracs[:KaiserCriterion])*100.0,2))+'% (based upon mean centered auto scaled input data) ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    LogString = '++ Performing Horn\'s parallel analysis to determine the number of factors needed ...'
    print LogString
    Log.Write(LogString+'\n')
    NPCs = HornsParallalAnalysis(PhenotypeArrayDict,
                                 Log)
    LogString  = '  ** Number of retained PCs from parallel analysis = '+str(NPCs)+' ...\n'
    LogString += '  ** The percentage explained variance for this number of PCs is '+str(round(scipy.sum(MyPCA.fracs[:NPCs])*100.0,2))+'% (based upon mean centered auto scaled input data) ...'
    print LogString
    Log.Write(LogString+'\n')
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    LogString = '++ 2D-plotting scores and loadings for all combinations of the retained number of PCs...'
    print LogString
    Log.Write(LogString+'\n')
    PlotScoresAndLoadings(MyData,
                          13,
#                          NPCs,
                          MyPCA,
                          MtbClassColors,
                          MtbClasses,
                          MtbName2ClassesDict,
                          MetaboliteList)
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