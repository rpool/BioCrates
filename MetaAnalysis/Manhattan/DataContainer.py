import re
import gzip
import scipy
import pylab

import Logger

class DataContainer:
    def __init__(self):
        self.DataName  = None
        self.DataArray = None
        return

    def SetDataName(self,
                    Name=str):
        self.DataName = Name
        return

    def GetDataName(self):
        return self.DataName

    def InitDataArray(self):
        self.DataArray = []
        return

    def AppendToArray(self,
                      Entry=str):
        self.DataArray.append(Entry)
        return

    def GetDataArray(self):
        return self.DataArray


class DataContainers:
    def __init__(self):
        self.Name2Column    = {}
        self.Column2Name    = {}
        self.DataContainers = {}
        self.Label          = ''
        self.Color          = None
        return

    def ParseGWAOutput(self,
                       FileName=str,
                       Log=Logger):
        LogString = '**** Parsing GWA output file \"'+FileName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')

        fr = None
        if(re.search('.gz',FileName)):
            fr = gzip.open(FileName,'r')
        else:
            fr = open(FileName,'r')

        Header      = fr.readline()
        CountColumn = 0
        for Name in Header.strip().split():
            self.DataContainers[Name]     = DataContainer()
            self.Column2Name[CountColumn] = Name
            self.Name2Column[Name]        = CountColumn
            self.DataContainers[Name].SetDataName(Name)
            self.DataContainers[Name].InitDataArray()
            CountColumn                  += 1
        for Line in fr:
            LSplit = Line.strip().split()
            for i in range(len(LSplit)):
                Entry = LSplit[i]
                Name  = self.Column2Name[i]
                self.DataContainers[Name].AppendToArray(Entry)

        if(type(fr)==file):
            fr.close()

        return

    def PylabUpdateParams(self):
        figwidth_pt   = 948.0 # pt (from revtex \showthe\columnwidth)
#        figwidth_pt   = 246.0 # pt (from revtex \showthe\columnwidth)
        inches_per_pt = 1.0/72.27
        figwidth      = figwidth_pt*inches_per_pt
        golden_mean   = (scipy.sqrt(5.0)-1.0)/2.0 # Aesthetic ratio
        figheight     = figwidth*golden_mean
        fig_size      = [figwidth,figheight]
        params        = {'backend': 'pdf',
                         'patch.antialiased': True,
                         'axes.labelsize': 10,
                         'axes.linewidth': 0.5,
                         'grid.color': '0.75',
                         'grid.linewidth': 0.25,
                         'grid.linestyle': ':',
                         'axes.axisbelow': False,
                         'text.fontsize': 10,
                         'legend.fontsize': 5,
                         'xtick.labelsize': 10,
                         'ytick.labelsize': 10,
                         'text.usetex': True,
                         'figure.figsize': fig_size}
        left   = 0.16
        bottom = 0.16
        width  = 0.86-left
        height = 0.95-bottom

        return params,\
               [left, bottom, width, height]

    def PlotManhattan(self,
                      xname=str,
                      yname=str,
                      Log=Logger):

        LogString = '**** Generating Manhattan plot ...'
        print LogString
        Log.Write(LogString+'\n')

        XName = ''
        YName = ''
        for Key in self.DataContainers.iterkeys():
            if(re.search(xname,Key)):
                XName = Key
            if(re.search(yname,Key)):
                YName = Key

        X = []
        for Entry in self.DataContainers[XName].GetDataArray():
            X.append(int(Entry))
        X = scipy.array(X)
        Y = []
        for Entry in self.DataContainers[YName].GetDataArray():
            if(Entry=='-1'):
                Y.append(float(1.0))
            else:
                Y.append(float(Entry))
        Y = -scipy.log10(scipy.array(Y))

        PylabParameters,\
        Rectangle         = self.PylabUpdateParams()
        pylab.rcParams.update(PylabParameters)
        PylabFigure = pylab.figure()
        PylabFigure.clf()
        PylabAxis = PylabFigure.add_axes(Rectangle)
        PylabAxis.scatter(x=X,
                          y=Y)
        XSign = scipy.array(PylabAxis.get_xlim())
        YSign = -scipy.log10(scipy.array([5.0e-8,5.0e-8]))
        PylabAxis.plot(XSign,
                       YSign,
                       linestyle='--',
                       color='grey',
                       linewidth=1.0)
        XSugg = scipy.array(PylabAxis.get_xlim())
        YSugg = -scipy.log10(scipy.array([1.0e-5,1.0e-5]))
        PylabAxis.plot(XSugg,
                       YSugg,
                       linestyle=':',
                       color='grey',
                       linewidth=1.0)
        PylabAxis.set_ylim([0.0,PylabAxis.get_ylim()[1]])
        PylabFigure.savefig('Manhattan.png')
        PylabAxis.clear()
        pylab.close(PylabFigure)
        del PylabFigure
        del PylabAxis

        return

class ListDataContainers:
    def __init__(self):
        self.List              = []
        self.PhenotypeName     = ''
        self.OffsetBetweenChrs = 0
        return

    def SetPhenotypeName(self,
                         Name=str):
        self.PhenotypeName = Name
        return

    def GetPhenotypeName(self):
        return self.PhenotypeName


    def PylabUpdateParams(self):
        figwidth_pt   = 1422.0 # pt (from revtex \showthe\columnwidth)
#        figwidth_pt   = 246.0 # pt (from revtex \showthe\columnwidth)
        inches_per_pt = 1.0/72.27
        figwidth      = figwidth_pt*inches_per_pt
        golden_mean   = (scipy.sqrt(5.0)-1.0)/2.0 # Aesthetic ratio
        figheight     = figwidth*golden_mean
        fig_size      = [figwidth,figheight]
        params        = {'backend': 'pdf',
                         'patch.antialiased': True,
                         'axes.labelsize': 18,
                         'axes.linewidth': 0.5,
                         'grid.color': '0.75',
                         'grid.linewidth': 0.25,
                         'grid.linestyle': ':',
                         'axes.axisbelow': False,
                         'text.fontsize': 14,
                         'legend.fontsize': 14,
                         'xtick.labelsize': 14,
                         'ytick.labelsize': 14,
                         'text.usetex': True,
                         'figure.figsize': fig_size}
        left   = 0.06
        bottom = 0.10
        width  = 0.95-left
        height = 0.95-bottom

        return params,\
               [left, bottom, width, height]

    def PlotManhattan(self,
                      xname=str,
                      yname=str,
                      Log=Logger):

        LogString = '**** Generating Manhattan plot ...'
        print LogString
        Log.Write(LogString+'\n')

        PylabParameters,\
        Rectangle         = self.PylabUpdateParams()
        pylab.rcParams.update(PylabParameters)
        PylabFigure = pylab.figure()
        PylabFigure.clf()
        PylabAxis   = PylabFigure.add_axes(Rectangle)
        XMax        = 0
        XTicks      = []
        XTickLabels = []
        for i in range(len(self.List)):
            DCs   = self.List[i]
            XName = ''
            YName = ''
            for Key in DCs.DataContainers.iterkeys():
                if(re.search(xname,Key)):
                    XName = Key
                if(re.search(yname,Key)):
                    YName = Key

            X = []
            Y = []
            for i in range(len(DCs.DataContainers[YName].GetDataArray())):
                YEntry = DCs.DataContainers[YName].GetDataArray()[i]
                XEntry = DCs.DataContainers[XName].GetDataArray()[i]
                if(YEntry!='-1'):
                    Y.append(float(YEntry))
                    X.append(int(XEntry))
                else:
                    Y.append(1.0)
                    X.append(int(XEntry))
            XTicks.append(0.5*(min(X)+max(X))+XMax)
            XTickLabels.append(r'${\rm '+DCs.Label+'}$')
            X    = scipy.array(X) + XMax
            XMax = X.max()+self.OffsetBetweenChrs
            Y    = -scipy.log10(scipy.array(Y))

            YInsign = Y < -scipy.log10(1.0e-6)
            YSign   = Y > -scipy.log10(5.0e-8)
            YSugg   = Y >= -scipy.log10(1.0e-6)
            YSugg  *= Y <= -scipy.log10(5.0e-8)

            YY = scipy.compress(YInsign,Y)
            if(len(YY)>0):
                PylabAxis.scatter(x=scipy.compress(YInsign,X),
                                  y=YY,
                                  color=DCs.Color,
                                  s=0.5)
            YY = scipy.compress(YSugg,Y)
            if(len(YY)>0):
                PylabAxis.scatter(x=scipy.compress(YSugg,X),
                                  y=YY,
                                  color=DCs.Color,
                                  s=5.0)
            YY = scipy.compress(YSign,Y)
            if(len(YY)>0):
                PylabAxis.scatter(x=scipy.compress(YSign,X),
                                  y=YY,
                                  color=DCs.Color,
                                  s=10.0)

        XSign = scipy.array(PylabAxis.get_xlim())
        YSign = -scipy.log10(scipy.array([5.0e-8,5.0e-8]))
        PylabAxis.plot(XSign,
                       YSign,
                       linestyle='--',
                       color='grey',
                       label=r'${\rm '+self.PhenotypeName+'}$',
                       linewidth=1.25)
        XSugg = scipy.array(PylabAxis.get_xlim())
        YSugg = -scipy.log10(scipy.array([1.0e-6,1.0e-6]))
        PylabAxis.plot(XSugg,
                       YSugg,
                       linestyle=':',
                       color='grey',
                       linewidth=1.25)
        PylabAxis.set_ylim([0.0,PylabAxis.get_ylim()[1]])
#        PylabAxis.set_xlim([-5e7,PylabAxis.get_xlim()[1]])
        PylabAxis.set_xlim([0,XMax])
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
        PylabAxis.xaxis.set_ticks(XTicks)
        PylabAxis.xaxis.set_ticklabels(XTickLabels)
        for Label in PylabAxis.xaxis.get_ticklabels():
            Label.set_rotation(90)
        PylabAxis.set_ylabel(r'$-{\rm log}_{10}(p-{\rm value})$')
        PylabFigure.savefig('Manhattan_'+self.PhenotypeName+'.png')
        PylabAxis.clear()
        pylab.close(PylabFigure)
        del PylabFigure
        del PylabAxis

        return