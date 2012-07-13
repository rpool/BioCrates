import os
import pylab
import scipy
import scipy.stats
import scipy.optimize

def LinearModel(X=scipy.array,
                Parameters=[]):
#   The linear regression function with ofsset 0.
    B = Parameters[0]
    Y = B*X
    return Y

def Residuals(Parameters=[],
              X=scipy.array,
              Y=scipy.array):
#   The residuals function to be called by scipy.optimize.leastsq().
    Err = Y - LinearModel(X,
                          Parameters)
    return Err

def GetPPointsArray(ArraySize=int):
    #===========================================================================
    # Copied fromm R-ppoints
    #===========================================================================
    PPointsArray   = [] # probability points array
    ThreeOverEight = 3.0/8.0
    Half           = 0.5
    N              = None
    Factor         = None
    if(ArraySize>1):
        N = ArraySize
    if(N<=10):
        Factor = ThreeOverEight
    else:
        Factor = Half
    OffSet       = Factor/float(N)
    PPointsArray = scipy.linspace(0.0+OffSet,1.0-OffSet,N)

    return PPointsArray

def LambdaEstimate(Array=scipy.array,
                   Filter=True):
    #===========================================================================
    # copied from R-GenABEL.lamdaest
    #===========================================================================
    Array      = Array.astype(float)
    Estimate   = None
    Ntp        = len(Array)
    QChi2Array = None

    if(Array.max()<=1.0):
#       Convert to quantile function of PValObsArray (df=1) if input are p-values.
        QChi2Array = scipy.stats.chi2.isf(Array,
                                          1)
    else:
        QChi2Array = Array

    if(Filter):
        FilterArray        = (QChi2Array>=1.0e-8)
        QChi2FilteredArray = scipy.compress(FilterArray,QChi2Array)
    else:
        QChi2FilteredArray = QChi2Array
    QChi2FilteredArray.sort()

    PPointsArray = GetPPointsArray(len(QChi2FilteredArray))
    PPointsArray = scipy.sort(scipy.stats.chi2.ppf(PPointsArray,1)) # df=1

    FilterArray          = (PPointsArray!=0.0)
    FilterArray         *= (QChi2FilteredArray!=0.0)
    PPointsArray         = scipy.compress(FilterArray,PPointsArray)
    QChi2FilteredArray   = scipy.compress(FilterArray,QChi2FilteredArray)

#   Fit PPointsArray,QChi2FilteredArray to the linear model.
    P0           = [1.0]
    PBest        = scipy.optimize.leastsq(Residuals,
                                          P0,
                                          args=(PPointsArray,QChi2FilteredArray),
                                          full_output=1,
                                          maxfev=100)
    Estimate = None
    if(type(PBest[0])==scipy.float64):
        Estimate = PBest[0]
    else:
        Estimate     = PBest[0][0]
#   Error estimation of parameter.
    Chi2 = scipy.power(PBest[2]['fvec'],2.0).sum()
    Dof  = len(QChi2FilteredArray)-len(P0)-1
    SE   = scipy.real(scipy.sqrt(PBest[1][0,0])*scipy.sqrt(Chi2/float(Dof)))

    return Estimate,SE

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
                   'grid.color': '0.75',
                   'grid.linewidth': 0.25,
                   'grid.linestyle': ':',
                   'axes.axisbelow': False,
                   'text.fontsize': 8,
                   'legend.fontsize': 4,
                   'xtick.labelsize': 8,
                   'ytick.labelsize': 8,
                   'text.usetex': True,
                   'figure.figsize': FigSize}
    Left   = 0.16
    Bottom = 0.16
    Width  = 0.86 - Left
    Height = 0.95 - Bottom

    return Params,\
           [Left, Bottom, Width, Height]

def Scatter(X=scipy.array,
            Y=scipy.array,
            Min=0.0,
            Max=1.0,
            Legend='',
            Name=str,
            XLabel='',
            YLabel='',
            boDrawY_EQ_X=False):

    PylabParameters,\
    Rectangle         = PylabGetParams()
    Size              = 1.5
    LineWidth         = 0.5
    pylab.rcParams.update(PylabParameters)
    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis = PylabFigure.add_axes(Rectangle)
    PlotName  = Name
    PlotPath  = os.path.join(os.getcwd(),'Plots')
    if(not os.path.isdir(PlotPath)):
        os.mkdir(PlotPath)
    PlotName  = os.path.join(PlotPath,PlotName)
    PylabAxis.scatter(X,
                      Y,
                      color='blue',
                      s=Size,
                      linewidths=LineWidth,
                      facecolor='None',
                      label=Legend)
    PylabAxis.plot([Min,Max],
                   [Min,Max],
                   color='black',
                   linestyle='--',
                   linewidth=LineWidth,
                   label=r'$y=x$')
    PylabAxis.set_ylim([Min,Max])
    PylabAxis.set_xlim([Min,Max])
    PylabAxis.set_xlabel(XLabel)
    PylabAxis.set_ylabel(YLabel)
    PylabAxis.spines['right'].set_visible(False)
    PylabAxis.spines['top'].set_visible(False)
    PylabAxis.xaxis.set_ticks_position('bottom')
    PylabAxis.yaxis.set_ticks_position('left')
    Handles,Labels = PylabAxis.get_legend_handles_labels()
    PylabAxis.legend(Handles,
                     Labels,
                     fancybox=True,
                     shadow=True,
                     loc='lower right')
    PylabAxis.grid(True)
    PylabFigure.savefig(PlotName,dpi=600)

    pylab.close(PylabFigure)
    del PylabAxis
    del PylabFigure

    return

def PlotQQFilteredOnScore(MtbName=str,
                          LPValExpArray=scipy.array,
                          LPValObsArray=scipy.array,
                          PValObsArray=scipy.array,
                          EMACArray=scipy.array,
                          EMACLevels=[],
                          Colors=[]):
    PylabParameters,\
    Rectangle         = PylabGetParams()
    Size              = 1.5
    LineWidth         = 0.5
    pylab.rcParams.update(PylabParameters)
    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis = PylabFigure.add_axes(Rectangle)
    PlotName  = 'QQPValueModeFilteredOnScore_'+MtbName+'.png'
    PlotPath  = os.path.join(os.getcwd(),'Plots')
    if(not os.path.isdir(PlotPath)):
        os.mkdir(PlotPath)
    PlotName = os.path.join(PlotPath,PlotName)
#   Initialize the Lamba arrays and start filling them (Lambda = genomic inflation factor)
    Lambdas   = []
    SEsLambda = []
    LambdaEst,\
    SELambdaEst = LambdaEstimate(Array=PValObsArray,Filter=True)
    Lambdas.append(LambdaEst)
    SEsLambda.append(SELambdaEst)
    #   Plot the unfiltered QQ plot.
    PylabAxis.scatter(scipy.sort(LPValExpArray),
                      scipy.sort(LPValObsArray),
                      color=Colors[0],
                      s=Size,
                      linewidths=LineWidth,
                      facecolor='None',
                      label=r'${\rm ~all~SNPs}'+r'~(\lambda='+str(round(LambdaEst,2))+r')$')
    for i in range(1,len(EMACLevels)):
#       Plot the Score filtered QQ plots
        Color       = Colors[i]
        LabelString = None

        FilterArray  = (EMACArray >= EMACLevels[i-1])
        FilterArray *= (EMACArray  < EMACLevels[i])
        ObsLP = scipy.sort(scipy.compress(FilterArray,LPValObsArray))
        ObsP  = scipy.sort(scipy.compress(FilterArray,PValObsArray))
#       Append estimates to Lambda arrays
        LambdaEst,\
        SELambdaEst = LambdaEstimate(Array=ObsP,Filter=True)
        Lambdas.append(LambdaEst)
        SEsLambda.append(SELambdaEst)
        ExpLP = scipy.sort(scipy.compress(FilterArray,LPValExpArray))

        LabelString  = str(EMACLevels[i-1])+r'\leq {\rm EMAC} <'+str(EMACLevels[i])
        LabelString += r'~(\lambda = '+str(round(LambdaEst,2))+')'
        PylabAxis.scatter(ExpLP,
                          ObsLP,
                          color=Color,
                          s=Size,
                          linewidths=LineWidth,
                          facecolor='None',
                          label=r'${\rm ~'+LabelString+'}$')
    MaxLPVal    = LPValObsArray.max()
    MaxLPValExp = LPValExpArray.max()
    Max         = max(MaxLPVal,MaxLPValExp)+2.5
    PylabAxis.plot([0.0,Max],
                   [0.0,Max],
                   color='black',
                   linestyle='--',
                   linewidth=LineWidth)
    PylabAxis.set_ylim([0.0,Max])
    PylabAxis.set_xlim([0.0,Max])
    PylabAxis.set_xlabel(r'$-\log_{10}{(p-{\rm value})} {\rm ~(expected)~[-]}$')
    PylabAxis.set_ylabel(r'$-\log_{10}{(p-{\rm value})} {\rm ~(observed)~[-]}$')
    PylabAxis.spines['right'].set_visible(False)
    PylabAxis.spines['top'].set_visible(False)
    PylabAxis.xaxis.set_ticks_position('bottom')
    PylabAxis.yaxis.set_ticks_position('left')
    Handles,Labels = PylabAxis.get_legend_handles_labels()
    PylabAxis.legend(Handles,
                     Labels,
                     fancybox=True,
                     shadow=True,
                     loc='lower right')
    PylabAxis.grid(True)
    PylabFigure.savefig(PlotName,dpi=600)

    pylab.close(PylabFigure)
    del PylabAxis
    del PylabFigure

    return Lambdas[0],\
           SEsLambda[0]