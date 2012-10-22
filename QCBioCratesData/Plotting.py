import scipy.optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import re

import BioCratesAnalyticalRanges
import Logger

def Gauss(x,
          p):
    A, mu, sigma = p

    return A*scipy.exp(-(x-mu)**2.0/(2.0*sigma**2.0))

def Residuals(p,
              y,
              x):
    Err = y - Gauss(x,p)

    return Err

def PlotDistribution(Ax=pyplot.Axes,
                     DataArray=[],
                     LODValue=float,
                     LLOQValue=float,
                     ULOQValue=float,
                     MetaboliteName=str,
                     MetaboliteClass=str,
                     XLabel=str,
                     YLabel=str,
                     NBins=100,
                     PlotName='Fig.pdf',
                     Color=str,
                     boFit=False,
                     Label=None,
                     boPlotAnalRanges=False):
    DataArray  = scipy.real(scipy.array(DataArray))
    Mu         = DataArray.mean()
    Sig        = DataArray.std()
    NormFactor = 1.0
    POfC,C     = scipy.histogram(a=DataArray,
                                 bins=NBins,
                                 normed=True)
    POfC      = POfC/POfC.sum()
    CCentered = 0.5*(C[1:]+C[:-1])

    Ax.plot(CCentered,
            POfC,
            color=Color,
            linestyle='-',
            linewidth=0.25,
            label=Label)
    if(boFit):
        P0              = [NormFactor,Mu,Sig] # initial guesses
        PBest           = scipy.optimize.leastsq(Residuals,P0,args=(POfC,CCentered),full_output=1,maxfev=100)
        BestParams      = PBest[0]
        Sig             = BestParams[2]
        Mu              = BestParams[1]
        NormFactor      = BestParams[0]
#        cov_x           = pbest[1]
#        print 'Best fit parameters for Bolzmann probability distribution ',bestparams
#        print cov_x
        step    = 1.0e-3
#        XRange  = scipy.arange(DataArray.min(),DataArray.max(),step)
        XRange  = scipy.arange(max(0.0,Mu-7.5*Sig),Mu+7.5*Sig,step)
        DataFit = Gauss(XRange,BestParams)
        Ax.plot(XRange,DataFit,linewidth=0.5,color=Color,ls=':')

    if(boPlotAnalRanges):
        Y = scipy.array([0.0,POfC.max()])
        if(LODValue!=0.0):
            X = scipy.array([LODValue,LODValue])
            Ax.plot(X,
                    Y,
                    color='grey',
                    linestyle='-',
                    linewidth=0.5,
                    label=r'${\rm LOD}$')
        if(LLOQValue!=0.0):
            X = scipy.array([LLOQValue,LLOQValue])
            Ax.plot(X,
                    Y,
                    color='grey',
                    linestyle='--',
                    linewidth=0.5,
                    label=r'${\rm LLOQ}$')
#        if(ULOQValue!=0.0):
#            X = scipy.array([ULOQValue,ULOQValue])
#            Ax.plot(X,
#                    Y,
#                    color='grey',
#                    linestyle=':',
#                    linewidth=0.5,
#                    label=r'${\rm ULOQ}$')
    Ax.set_xlabel(r'$'+XLabel+'$')
    Ax.set_ylabel(r'$'+YLabel+'$')
    Ax.spines['right'].set_visible(False)
    Ax.spines['top'].set_visible(False)
    Ax.xaxis.set_ticks_position('bottom')
    Ax.yaxis.set_ticks_position('left')
    Handles,Labels = Ax.get_legend_handles_labels()
    Ax.legend(Handles,
              Labels,
              fancybox=True,
              shadow=True,
              loc="best")
    Ax.grid(True)

    return

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

def PlotDistributions(DataDict=dict,
                      AnalytRanges=BioCratesAnalyticalRanges,
                      PlotBaseName=str,
                      Extension='.pdf',
                      Resolution=400,
                      Log=Logger):
    LogString = '** Plotting distributions of measured metabolite concentrations ...'
    print LogString
    Log.Write(LogString+'\n')
    Counter    = 0
    PlotNames  = []
    PylabParameters,\
    Rectangle         = PylabUpdateParams()
    pyplot.rcParams.update(PylabParameters)
    for Value in DataDict.itervalues():
        if(Value.GetMetaboliteName()):
            Counter += 1
            LogString = '** Plotting metabolite '+\
                        str(Counter)+\
                        ' with name / class: \"'+\
                        Value.GetMetaboliteName()+\
                        '\" / \"'+\
                        Value.GetMetaboliteClass()+\
                        '\" ...'
            print LogString
            Log.Write(LogString+'\n')

            PylabFigure = pyplot.figure()
            PylabFigure.clf()
            PylabAxis   = PylabFigure.add_axes(Rectangle)

            PlotName  = re.sub(Extension,str(Counter)+Extension,PlotBaseName)
            PlotNames.append(PlotName)
            NBins = 30

            DataArray  = scipy.real(scipy.array(Value.GetDataArray()))
            POfC,C     = scipy.histogram(a=DataArray,
                                         bins=NBins,
                                         normed=True)
            POfC      = POfC/POfC.sum()
            CCentered = 0.5*(C[1:]+C[:-1])
            MinWidth  = 1.0e500
            for i in range(len(CCentered)-1):
                MinWidth = min(MinWidth,abs(CCentered[i+1]-CCentered[i]))
            CLeft = CCentered - 0.5*MinWidth

            LODValue = AnalytRanges.LODValueDict[str(Value.GetMetaboliteName())]
            Y = scipy.array([0.0,POfC.max()])
            if(LODValue!=None):
                X = scipy.array([LODValue,LODValue])
                PylabAxis.plot(X,
                               Y,
                               color='grey',
                               linestyle='-',
                               linewidth=0.5,
                               label=r'${\rm LOD}$')
            LLOQValue = AnalytRanges.LLOQValueDict[str(Value.GetMetaboliteName())]
            if(LLOQValue!=None):
                X = scipy.array([LLOQValue,LLOQValue])
                PylabAxis.plot(X,
                               Y,
                               color='grey',
                               linestyle='--',
                               linewidth=0.5,
                               label=r'${\rm LLOQ}$')

            PylabAxis.bar(left=CLeft,
                          height=POfC,
                          width=MinWidth,
                          bottom=0.0,
                          color='black',
                          linewidth=0.25,
                          alpha=0.15,
                          label=r'${\tt '+Value.GetMetaboliteConventionName()+r'}{\rm ~(raw)}$')

            Mu              = DataArray.mean()
            Sig             = DataArray.std()
            LogString       = '** Raw data mean / std        : '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
            NormFactor      = 1.0
            P0              = [NormFactor,Mu,Sig] # initial guesses
            PBest           = scipy.optimize.leastsq(Residuals,
                                                     P0,
                                                     args=(POfC,CCentered),
                                                     full_output=1,
                                                     maxfev=100)
            BestParams      = PBest[0]
            Sig             = BestParams[2]
            Mu              = BestParams[1]
            NormFactor      = BestParams[0]
            SSErr           = (PBest[2]['fvec']**2).sum()
            SSTot           = ((POfC-POfC.mean())**2).sum()
            R2              = 1.0-(SSErr/SSTot)
            LogString      += '** Raw data fit mean / std    : '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
            LogString      += '** Raw data fit RSquared      : '+str(round(R2,3))
            print LogString
            Log.Write(LogString+'\n')

            NStep   = 1000
            Range   = (Mu+7.5*Sig)-max(0.0,Mu-7.5*Sig)
            Step    = Range/float(NStep)
            XRange  = scipy.arange(max(0.0,Mu-7.5*Sig),Mu+7.5*Sig,Step)
            DataFit = Gauss(XRange,BestParams)
            PylabAxis.plot(XRange,
                           DataFit,
                           linewidth=0.5,
                           color='black',
                           ls=':')

            QCedDataArray = Value.GetQCedDataArray()
            if(QCedDataArray!=None):
                DelList = []
                for i in range(len(QCedDataArray)):
                    Entry = QCedDataArray[i]
                    if(Entry==Value.GetMissingDataIdentifier()):
                        DelList.append(i)
                DelList.sort()
                DelList.reverse()
                for i in DelList:
                    del QCedDataArray[i]
                POfC,C     = scipy.histogram(a=scipy.array(QCedDataArray),
                                             bins=NBins,
                                             normed=True)
                POfC      = POfC/POfC.sum()
                CCentered = 0.5*(C[1:]+C[:-1])
                MinWidth  = 1.0e500
                for i in range(len(CCentered)-1):
                    MinWidth = min(MinWidth,abs(CCentered[i+1]-CCentered[i]))
                CLeft = CCentered - 0.5*MinWidth
                PylabAxis.bar(left=CLeft,
                              height=POfC,
                              width=MinWidth,
                              bottom=0.0,
                              color='blue',
                              linewidth=0.25,
                              alpha=0.15,
                              label=r'${\tt '+Value.GetMetaboliteConventionName()+r'}{\rm ~(QCed)}$')

                Mu              = scipy.array(QCedDataArray).mean()
                Sig             = scipy.array(QCedDataArray).std()
                NormFactor      = 1.0
                LogString       = '** QCed data mean / std       : '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
                P0              = [NormFactor,Mu,Sig] # initial guesses
                PBest           = scipy.optimize.leastsq(Residuals,
                                                         P0,
                                                         args=(POfC,CCentered),
                                                         full_output=1,
                                                         maxfev=100)
                BestParams      = PBest[0]
                Sig             = BestParams[2]
                Mu              = BestParams[1]
                NormFactor      = BestParams[0]
                SSErr           = (PBest[2]['fvec']**2).sum()
                SSTot           = ((POfC-POfC.mean())**2).sum()
                R2              = 1.0-(SSErr/SSTot)
                LogString      += '** QCed data fit mean / std   : '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
                LogString      += '** QCed data fit RSquared     : '+str(round(R2,3))
                print LogString
                Log.Write(LogString+'\n')

                NStep   = 1000
                Range   = (Mu+7.5*Sig)-max(0.0,Mu-7.5*Sig)
                Step    = Range/float(NStep)
                XRange  = scipy.arange(max(0.0,Mu-7.5*Sig),Mu+7.5*Sig,Step)
                DataFit = Gauss(XRange,BestParams)
                PylabAxis.plot(XRange,
                               DataFit,
                               linewidth=0.5,
                               color='blue',
                               ls=':')

            ImputedDataArray = Value.GetImputedDataArray()
            if(ImputedDataArray!=None):
                POfC,C     = scipy.histogram(a=scipy.array(ImputedDataArray),
                                             bins=NBins,
                                             normed=True)
                POfC      = POfC/POfC.sum()
                CCentered = 0.5*(C[1:]+C[:-1])
                MinWidth  = 1.0e500
                for i in range(len(CCentered)-1):
                    MinWidth = min(MinWidth,abs(CCentered[i+1]-CCentered[i]))
                CLeft = CCentered - 0.5*MinWidth
                PylabAxis.bar(left=CCentered,
                              height=POfC,
                              width=MinWidth,
                              bottom=0.0,
                              color='red',
                              linewidth=0.25,
                              alpha=0.15,
                              label=r'${\tt '+Value.GetMetaboliteConventionName()+r'}{\rm ~(imputed)}$')

                Mu              = scipy.array(ImputedDataArray).mean()
                Sig             = scipy.array(ImputedDataArray).std()
                NormFactor      = 1.0
                LogString       = '** Imputed data mean / std    : '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
                P0              = [NormFactor,Mu,Sig] # initial guesses
                PBest           = scipy.optimize.leastsq(Residuals,
                                                         P0,
                                                         args=(POfC,CCentered),
                                                         full_output=1,
                                                         maxfev=100)
                BestParams      = PBest[0]
                Sig             = BestParams[2]
                Mu              = BestParams[1]
                NormFactor      = BestParams[0]
                SSErr           = (PBest[2]['fvec']**2).sum()
                SSTot           = ((POfC-POfC.mean())**2).sum()
                R2              = 1.0-(SSErr/SSTot)
                LogString      += '** Imputed data fit mean / std: '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
                LogString      += '** Imputed data fit RSquared  : '+str(round(R2,3))
                print LogString
                Log.Write(LogString+'\n')

                NStep   = 1000
                Range   = (Mu+7.5*Sig)-max(0.0,Mu-7.5*Sig)
                Step    = Range/float(NStep)
                XRange  = scipy.arange(max(0.0,Mu-7.5*Sig),Mu+7.5*Sig,Step)
                DataFit = Gauss(XRange,BestParams)
                PylabAxis.plot(XRange,
                               DataFit,
                               linewidth=0.5,
                               color='red',
                               ls=':')

            XLabel = r'c~\rm [\mu M]'
            YLabel = r'P(c)~\rm [-]'
            PylabAxis.set_xlabel(r'$'+XLabel+'$')
            PylabAxis.set_ylabel(r'$'+YLabel+'$')
            PylabAxis.spines['right'].set_visible(False)
            PylabAxis.spines['top'].set_visible(False)
            PylabAxis.xaxis.set_ticks_position('bottom')
            PylabAxis.yaxis.set_ticks_position('left')
            Handles,Labels = PylabAxis.get_legend_handles_labels()
            PylabAxis.legend(Handles,
                             Labels,
                             fancybox=True,
                             shadow=True,
                             loc="best")
            PylabAxis.grid(True)
            PylabFigure.savefig(PlotName,
                                dpi=Resolution)
            PylabAxis.clear()
            pyplot.close(PylabFigure)
            del PylabFigure
            del PylabAxis

    return PlotNames

def PlotLnSpaceDistributions(DataDict=dict,
                             AnalytRanges=BioCratesAnalyticalRanges,
                             PlotBaseName=str,
                             Extension='.pdf',
                             Resolution=400,
                             Log=Logger):
    LogString = '** Plotting distributions of ln-transformed measured metabolite concetrations ...'
    print LogString
    Log.Write(LogString+'\n')
    Counter    = 0
    PlotNames  = []
    PylabParameters,\
    Rectangle         = PylabUpdateParams()
    pyplot.rcParams.update(PylabParameters)
    for Value in DataDict.itervalues():
        if(Value.GetMetaboliteName()):
            Counter += 1
            LogString = '** Plotting metabolite '+\
                        str(Counter)+\
                        ' with name / class: \"'+\
                        Value.GetMetaboliteName()+\
                        '\" / \"'+\
                        Value.GetMetaboliteClass()+\
                        '\" ...'
            print LogString
            Log.Write(LogString+'\n')

            PylabFigure = pyplot.figure()
            PylabFigure.clf()
            PylabAxis   = PylabFigure.add_axes(Rectangle)

            PlotName  = re.sub(Extension,str(Counter)+Extension,PlotBaseName)
            PlotNames.append(PlotName)
            NBins = 30

            DataArray    = scipy.real(scipy.array(Value.GetDataArray()))
            LenDataArray = len(DataArray)
            FilterArray  = (DataArray!=0.0)
            DataArray    = scipy.compress(FilterArray,DataArray)
            DataArray    = scipy.log(DataArray)
            DLen         = LenDataArray - len(DataArray)
            if(DLen>0):
                LogString = '  !! Found '+str(DLen)+' 0.0 values in data array of metabolite \"'+Value.GetMetaboliteName()+'\" ...'
                print LogString
                Log.Write(LogString+'\n')
            POfC,C    = scipy.histogram(a=DataArray,
                                        bins=NBins,
                                        normed=True)
            POfC      = POfC/POfC.sum()
            CCentered = 0.5*(C[1:]+C[:-1])
            MinWidth  = 1.0e500
            for i in range(len(CCentered)-1):
                MinWidth = min(MinWidth,abs(CCentered[i+1]-CCentered[i]))
            CLeft = CCentered - 0.5*MinWidth

            LODValue = scipy.log(AnalytRanges.LODValueDict[str(Value.GetMetaboliteName())])
            Y = scipy.array([0.0,POfC.max()])
            if(LODValue!=None):
                X = scipy.array([LODValue,LODValue])
                PylabAxis.plot(X,
                               Y,
                               color='grey',
                               linestyle='-',
                               linewidth=0.5,
                               label=r'${\rm LOD}$')
            LLOQValue = scipy.log(AnalytRanges.LLOQValueDict[str(Value.GetMetaboliteName())])
            if(LLOQValue!=None):
                X = scipy.array([LLOQValue,LLOQValue])
                PylabAxis.plot(X,
                               Y,
                               color='grey',
                               linestyle='--',
                               linewidth=0.5,
                               label=r'${\rm LLOQ}$')

            PylabAxis.bar(left=CLeft,
                          height=POfC,
                          width=MinWidth,
                          bottom=0.0,
                          color='black',
                          linewidth=0.25,
                          alpha=0.15,
                          label=r'${\tt '+Value.GetMetaboliteConventionName()+r'}{\rm ~(raw)}$')

            Mu              = DataArray.mean()
            Sig             = DataArray.std()
            LogString       = '** Raw data mean / std        : '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
            NormFactor      = 1.0
            P0              = [NormFactor,Mu,Sig] # initial guesses
            PBest           = scipy.optimize.leastsq(Residuals,
                                                     P0,
                                                     args=(POfC,CCentered),
                                                     full_output=1,
                                                     maxfev=100)
            BestParams      = PBest[0]
            Sig             = BestParams[2]
            Mu              = BestParams[1]
            NormFactor      = BestParams[0]
            SSErr           = (PBest[2]['fvec']**2).sum()
            SSTot           = ((POfC-POfC.mean())**2).sum()
            R2              = 1.0-(SSErr/SSTot)
            LogString      += '** Raw data fit mean / std    : '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
            LogString      += '** Raw data fit RSquared      : '+str(round(R2,3))
            print LogString
            Log.Write(LogString+'\n')

            NStep   = 1000
            Range   = (Mu+7.5*Sig)-(Mu-7.5*Sig)
            Step    = Range/float(NStep)
            XRange  = scipy.arange(Mu-7.5*Sig,Mu+7.5*Sig,Step)
            DataFit = Gauss(XRange,BestParams)
            PylabAxis.plot(XRange,
                           DataFit,
                           linewidth=0.5,
                           color='black',
                           ls=':')

            QCedDataArray = Value.GetQCedDataArray()
            if(QCedDataArray!=None):
                DelList = []
                for i in range(len(QCedDataArray)):
                    Entry = QCedDataArray[i]
                    if((Entry==Value.GetMissingDataIdentifier()) or
                       (Entry==0.0)):
                        DelList.append(i)
                DelList.sort()
                DelList.reverse()
                for i in DelList:
                    del QCedDataArray[i]
                POfC,C      = scipy.histogram(a=scipy.log(scipy.array(QCedDataArray)),
                                              bins=NBins,
                                              normed=True)
                POfC      = POfC/POfC.sum()
                CCentered = 0.5*(C[1:]+C[:-1])
                MinWidth  = 1.0e500
                for i in range(len(CCentered)-1):
                    MinWidth = min(MinWidth,abs(CCentered[i+1]-CCentered[i]))
                CLeft = CCentered - 0.5*MinWidth
                PylabAxis.bar(left=CLeft,
                              height=POfC,
                              width=MinWidth,
                              bottom=0.0,
                              color='blue',
                              linewidth=0.25,
                              alpha=0.15,
                              label=r'${\tt '+Value.GetMetaboliteConventionName()+r'}{\rm ~(QCed)}$')

                Mu              = scipy.log(scipy.array(QCedDataArray)).mean()
                Sig             = scipy.log(scipy.array(QCedDataArray)).std()
                NormFactor      = 1.0
                LogString       = '** QCed data mean / std       : '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
                P0              = [NormFactor,Mu,Sig] # initial guesses
                PBest           = scipy.optimize.leastsq(Residuals,
                                                         P0,
                                                         args=(POfC,CCentered),
                                                         full_output=1,
                                                         maxfev=100)
                BestParams      = PBest[0]
                Sig             = BestParams[2]
                Mu              = BestParams[1]
                NormFactor      = BestParams[0]
                SSErr           = (PBest[2]['fvec']**2).sum()
                SSTot           = ((POfC-POfC.mean())**2).sum()
                R2              = 1.0-(SSErr/SSTot)
                LogString      += '** QCed data fit mean / std   : '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
                LogString      += '** QCed data fit RSquared     : '+str(round(R2,3))
                print LogString
                Log.Write(LogString+'\n')

                NStep   = 1000
                Range   = (Mu+7.5*Sig)-(Mu-7.5*Sig)
                Step    = Range/float(NStep)
                XRange  = scipy.arange(Mu-7.5*Sig,Mu+7.5*Sig,Step)
                DataFit = Gauss(XRange,BestParams)
                PylabAxis.plot(XRange,
                               DataFit,
                               linewidth=0.5,
                               color='blue',
                               ls=':')

            ImputedDataArray = Value.GetImputedDataArray()
            if(ImputedDataArray!=None):
                ImputedDataArray = scipy.array(ImputedDataArray)
                FilterArray      = (ImputedDataArray!=0.0)
                ImputedDataArray = scipy.compress(FilterArray,ImputedDataArray)

                POfC,C     = scipy.histogram(a=scipy.log(ImputedDataArray),
                                             bins=NBins,
                                             normed=True)
                POfC      = POfC/POfC.sum()
                CCentered = 0.5*(C[1:]+C[:-1])
                MinWidth  = 1.0e500
                for i in range(len(CCentered)-1):
                    MinWidth = min(MinWidth,abs(CCentered[i+1]-CCentered[i]))
                CLeft = CCentered - 0.5*MinWidth
                PylabAxis.bar(left=CCentered,
                              height=POfC,
                              width=MinWidth,
                              bottom=0.0,
                              color='red',
                              linewidth=0.25,
                              alpha=0.15,
                              label=r'${\tt '+Value.GetMetaboliteConventionName()+r'}{\rm ~(imputed)}$')

                Mu              = scipy.log(ImputedDataArray).mean()
                Sig             = scipy.log(ImputedDataArray).std()
                NormFactor      = 1.0
                LogString       = '** Imputed data mean / std    : '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
                P0              = [NormFactor,Mu,Sig] # initial guesses
                PBest           = scipy.optimize.leastsq(Residuals,
                                                         P0,
                                                         args=(POfC,CCentered),
                                                         full_output=1,
                                                         maxfev=100)
                BestParams      = PBest[0]
                Sig             = BestParams[2]
                Mu              = BestParams[1]
                NormFactor      = BestParams[0]
                SSErr           = (PBest[2]['fvec']**2).sum()
                SSTot           = ((POfC-POfC.mean())**2).sum()
                R2              = 1.0-(SSErr/SSTot)
                LogString      += '** Imputed data fit mean / std: '+str(round(Mu,3))+' / '+str(round(Sig,3))+'\n'
                LogString      += '** Imputed data fit RSquared  : '+str(round(R2,3))
                print LogString
                Log.Write(LogString+'\n')

                NStep   = 1000
                Range   = (Mu+7.5*Sig)-(Mu-7.5*Sig)
                Step    = Range/float(NStep)
                XRange  = scipy.arange(Mu-7.5*Sig,Mu+7.5*Sig,Step)
                DataFit = Gauss(XRange,BestParams)
                PylabAxis.plot(XRange,
                               DataFit,
                               linewidth=0.5,
                               color='red',
                               ls=':')

            XLabel = r'\ln{[c]}~\rm [-]'
            YLabel = r'P(c)~\rm [-]'
            PylabAxis.set_xlabel(r'$'+XLabel+'$')
            PylabAxis.set_ylabel(r'$'+YLabel+'$')
            PylabAxis.spines['right'].set_visible(False)
            PylabAxis.spines['top'].set_visible(False)
            PylabAxis.xaxis.set_ticks_position('bottom')
            PylabAxis.yaxis.set_ticks_position('left')
            Handles,Labels = PylabAxis.get_legend_handles_labels()
            PylabAxis.legend(Handles,
                             Labels,
                             fancybox=True,
                             shadow=True,
                             loc="best")
            PylabAxis.grid(True)
            PylabFigure.savefig(PlotName,
                                dpi=Resolution)
            PylabAxis.clear()
            pyplot.close(PylabFigure)
            del PylabFigure
            del PylabAxis

    return PlotNames