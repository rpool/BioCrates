#! /usr/bin/env python
import os
import sys
import re
import scipy
import pylab
import matplotlib

import ArgumentParser
import Logger
import fnmatch

iLARGE = 1e500
fLARGE = 1.0e500

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
    if(Arguments.YProperty=='PHE'):
        XPath               = os.path.join(Arguments.PositionPath,'SNPInfo.txt.gz')
        SNPInfoDecompressed = os.path.join(TmpPath,'SNPInfo.txt')
        os.system('pigz -d -c -k '+XPath+' > '+SNPInfoDecompressed)
        FH          = open(SNPInfoDecompressed,'r')
        HeaderList  = FH.readline().strip().split()
        FH.close()
        SNPIDArray  = None
        SNPIDList   = None
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
            XLeft.append(float(TmpTmpXPosArray.min())+XMax)
            XRight.append(float(TmpTmpXPosArray.max())+XMax)
            XTicks.append(0.5*(XLeft[-1]+XRight[-1]))
            XTickLabels.append(r'${\rm CHR'+str(Chr)+r'}$')

            TmpSNPIDArray = scipy.append(TmpSNPIDArray,TmpTmpSNPIDArray)
            TmpChrArray   = scipy.append(TmpChrArray,TmpTmpChrArray)
            TmpPosArray   = scipy.append(TmpPosArray,TmpTmpPosArray)
            TmpXPosArray  = scipy.append(TmpXPosArray,TmpTmpXPosArray)
            XMax          = TmpTmpXPosArray.max()
            del TmpTmpSNPIDArray
            del TmpTmpChrArray
            del TmpTmpPosArray
            del TmpTmpXPosArray

        SNPIDArray = scipy.array(TmpSNPIDArray)
        SNPIDList  = SNPIDArray.tolist()
        ChrArray   = scipy.array(TmpChrArray)
        PosArray   = scipy.array(TmpPosArray)
        XPosArray  = scipy.array(TmpXPosArray)
        del TmpSNPIDArray
        del TmpChrArray
        del TmpPosArray
        del TmpXPosArray

        for i in range(len(SNPIDArray)):
            Entry                      = SNPIDArray[i]
            SNPInfoDict[Entry]         = {}
            SNPInfoDict[Entry]['chr']  = ChrArray[i]
            SNPInfoDict[Entry]['pos']  = PosArray[i]
            SNPInfoDict[Entry]['xpos'] = XPosArray[i]
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

#    if(Arguments.YProperty=='PHE'):
#        XPath       = os.path.join(Arguments.SnpTestPath,Arguments.XProperty)
#        XListDir    = os.listdir(XPath)
#        MinXSpacing = iLARGE
#        XMin        = None
#        XMax        = 0
#        X           = []
#        XLeft       = []
#        XRight      = []
#        XTicks      = []
#        XTickLabels = []
#        for p in range(1): # if XProperty=='pos', all phenotypes have the same pos file content.
#            P = 'PHE'+str(p+1)+'_'
#            for c in range(Arguments.NChr):
#                C = 'CHR'+str(c+1)+'_'
#                for File in XListDir:
#                    if(re.search(C+P,File)):
#                        fr   = open(os.path.join(XPath,File),'r')
#                        FMem = []
#                        for Line in fr.readlines():
#                            FMem.append(int(Line.strip().split()[0]))
#                        fr.close()
#                        for i in range(1,len(FMem)):
#                            MinXSpacing = min(MinXSpacing,FMem[i]-FMem[i-1])
##                        XMin = min(FMem)+XMax
#                        X.extend(list(scipy.array(FMem)+XMax))
#                        if(p==0):
#                            XLeft.append(float(min(FMem))+XMax)
#                            XRight.append(float(max(FMem))+XMax)
#                            XTicks.append(0.5*(float(min(FMem))+float(max(FMem)))+XMax)
#                            XTickLabels.append(r'${\rm '+re.sub('_','',C)+'}$')
#                        XMax = max(X)
#                        del FMem
#        X     = scipy.array(X)
#        XXMin = X.min()
#        XXMax = X.max()

        LogString = '**** Parsed x-axis properties ...'
        print LogString
        Log.Write(LogString+'\n')

        YMin        = 1
        YMax        = Arguments.NPhe
        MinYSpacing = 1
        Y           = []
        for y in range(YMax):
            Y.append(y+1)
        LogString = '**** Processed y-axis properties ...'
        print LogString
        Log.Write(LogString+'\n')

        MNumber = None
        MName   = None
        MClass  = None
        if(Arguments.MetabolitClassesFileName!=''):
            MNumber = []
            MName   = []
            MClass  = []
            fr      = open(Arguments.MetabolitClassesFileName,'r')
            for Line in fr.readlines():
                if(Line[0]=="#"):
                    continue
                LStrp = Line.strip().split()
                MNumber.append(LStrp[0])
                MName.append(LStrp[1])
                MClass.append(LStrp[2])
            fr.close()
            Classes    = list(set(MClass))
            ClassDict  = {}
            ClassRange = {}
            for Entry in Classes:
                ClassDict[Entry] = []
                for i in range(len(MClass)):
                    if(Entry==MClass[i]):
                        ClassDict[Entry].append(int(MNumber[i]))
                ClassRange[Entry] = [min(ClassDict[Entry]),max(ClassDict[Entry])]

        ZMin        = 0.0
        ZMax        = -fLARGE
        ZPath       = os.path.join(Arguments.MAOutputPath,Arguments.ZProperty)
        ZListDir    = os.listdir(ZPath)
        if(Arguments.boReadZExtrFromFile):
            LogString = '**** Reading z-axis extrema from file \"'+os.path.join(ZPath,'Extrema.dat')+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            fr    = open(os.path.join(ZPath,'Extrema.dat'),'r')
            LSplt = fr.readline().strip().split()
            fr.close()
            ZMin = float(LSplt[0])
            ZMax = float(LSplt[1])
#        else:
#            for p in range(Arguments.NPhe): # Scan for the maximum value of Z for the colorbar
#                P = 'PHE'+str(p+1)+'_'
#                Z = []
#                for c in range(Arguments.NChr):
#                    C    = 'CHR'+str(c+1)+'_'
#                    File = os.path.join(ZPath,C+P+Arguments.ZProperty+'.dat')
#                    fr   = None
#                    if(os.path.basename(File) in ZListDir):
#                        fr = open(File,'r')
#                        for Line in fr:
#                            Value = Line.strip()
#                            if((Value!='-1') or
#                               (not re.search('nan',Value))):
#                                Z.append(float(Value))
#                            else:
#                                Z.append(1.0)
#                        fr.close()
#                Z    = scipy.array(Z)
#                ZMax = max(ZMax,scipy.real(-scipy.log10(Z)).max())
#                del Z
#                LogString = '**** Determining maximal z-axis value, now at '+re.sub('_','',P)+', ZMax = '+str(ZMax)+' ...'
#                print LogString
#                Log.Write(LogString+'\n')
#            LogString = '**** Writing z-axis extrema to file \"'+os.path.join(ZPath,'Extrema.dat')+'\" ...'
#            print LogString
#            Log.Write(LogString+'\n')
#            fw = open(os.path.join(ZPath,'Extrema.dat'),'w')
#            fw.write(str(ZMin)+' '+str(ZMax)+'\n')
#            fw.close()

        LogString = '**** Parsing in z-axis properties ...'
        print LogString
        Log.Write(LogString+'\n')
        ZSugg = {}
        ZSign = {}
        YSugg = {}
        YSign = {}
        XSugg = {}
        XSign = {}
        for p in range(Arguments.NPhe):
#        for p in range(0):
            P          = 'PHE'+str(p+1)+'_'
            PHE        = re.sub('_','',P)
            Mtb        = MName[MNumber.index(str(p+1))]
            FName      = os.path.join(Arguments.MAOutputPath,'MetaAnalysis_'+Mtb+'_1.tbl')
            if((not os.path.isfile(FName)) or
               (not os.path.islink(FName))):
                continue
            LogString  = '** Now at '+PHE+' (\"'+FName+'\") ...'
            print LogString
            Log.Write(LogString+'\n')
            FH     = open(FName,'r')
            Header = FH.readline().strip().split()
            FH.close()
            SNPIDCol = Header.index('MarkerName')
            PValCol  = Header.index('P-value')
            Arrays   = scipy.loadtxt(fname=FName,
                                     dtype=str,
                                     skiprows=1,
                                     usecols=[SNPIDCol,PValCol],
                                     unpack=True)
            RSIdArray  = Arrays[0]
            PValArray  = Arrays[1].astype(float)
            ZZ         = scipy.real(-scipy.log10(PValArray))
            TmpArray   = scipy.append(SNPIDArray,RSIdArray)
            TmpArray,\
            IndexArray = scipy.unique(ar=TmpArray,
                                      return_inverse=True)
            TmpArray   = scipy.append(RSIdArray,SNPIDArray)
            TmpArray,\
            tTmpArray  = scipy.unique(ar=TmpArray,
                                      return_inverse=True)
            IndexArray = scipy.append(IndexArray,tTmpArray)
            del TmpArray
            del tTmpArray
            IndexArray = scipy.unique(ar=IndexArray)
#            for Entry in RSIdArray:
#                IndexArray.append(SNPIDList.index(Entry))
#                print Entry
#            IndexArray = scipy.array(IndexArray)
            X          = XPosArray[IndexArray]
            ChromArray = ChrArray[IndexArray]
#            for c in range(Arguments.NChr):
#                C = 'CHR'+str(c+1)+'_'
#                for File in ZListDir:
#                    if(re.search(C+P,File)):
#                        fr   = open(os.path.join(ZPath,File),'r')
#                        FMem = []
#                        for Line in fr.readlines():
#                            z = Line.strip().split()[0]
#                            if(z!='-1'):
#                                FMem.append(float(z))
#                            else:
#                                FMem.append(1.0)
#                        fr.close()
#                        FMem  = scipy.real(-scipy.log10(scipy.array(FMem)))
#                        ZZ.extend(list(FMem))
#                        del FMem
            Sign  = ZZ  > -scipy.log10(5.0e-8)/float(Arguments.NPhe)
            Sugg  = ZZ >= -scipy.log10(1.0e-6)/float(Arguments.NPhe)
            Sugg *= ZZ <= -scipy.log10(5.0e-8)/float(Arguments.NPhe)
            ZSugg[PHE] = scipy.compress(Sugg,ZZ)
            ZSign[PHE] = scipy.compress(Sign,ZZ)
            YSugg[PHE] = scipy.ones(len(ZSugg[PHE]))*(p+1)
            YSign[PHE] = scipy.ones(len(ZSign[PHE]))*(p+1)
            XSugg[PHE] = scipy.compress(Sugg,X)
            XSign[PHE] = scipy.compress(Sign,X)
            del ZZ
            LogString  = '** Parsed \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')

        del X
        del Y

        LogString = '**** Generating 3D-Manhattan plot ...'
        print LogString
        Log.Write(LogString+'\n')

#       Lay-out stuff
        figwidth_pt   = 1422.0 # pt (from revtex \showthe\columnwidth)
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
        width  = 0.92-left
        height = 0.95-bottom

        pylab.rcParams.update(params)

        PylabFigure = pylab.figure()
        PylabFigure.clf()

        Rectangle   = [left, bottom, width, height]
        PylabAxis   = PylabFigure.add_axes(Rectangle)
        if(Arguments.XProperty=='pos'):
            PylabAxis.set_xlabel(r'${\rm position}$')
        if(Arguments.YProperty=='PHE'):
            PylabAxis.set_ylabel(r'${\rm metabolite}$')

        for p in range(Arguments.NPhe):
#        for p in range(0):
            P  = 'PHE'+str(p+1)
            if(len(ZSugg[P])>0):
                PylabAxis.scatter(x=XSugg[P],
                                  y=YSugg[P],
                                  color='black',
                                  s=ZSugg[P]/ZMax*100.0,
                                  marker='s',
                                  alpha=0.15,
                                  antialiased=True,
                                  edgecolors='none')
            if(len(ZSign[P])>0):
                PylabAxis.scatter(x=XSign[P],
                                  y=YSign[P],
                                  color='black',
                                  s=ZSign[P]/ZMax*100.0,
                                  marker='o',
                                  alpha=0.15,
                                  antialiased=True,
                                  edgecolors='none')
        PlotFile  = 'Manhattan3D.pdf'
        LogString = '**** Writing plot to \"'+PlotFile+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        XXRange  = float(XXMax)-float(XXMin)
        XXOffset = XXRange*0.005
        PylabAxis.set_xlim([float(XXMin)-XXOffset,float(XXMax)+XXOffset])
        PylabAxis.set_ylim([0,YMax+2])
        for Key,Value in ClassRange.iteritems():
            PylabAxis.plot(scipy.array([float(XXMin)-XXOffset,float(XXMax)+XXOffset]),
                           scipy.ones(2)*Value[0]-0.5,
                           lw=0.25,
                           color='black')
            PylabAxis.plot(scipy.array([float(XXMin)-XXOffset,float(XXMax)+XXOffset]),
                           scipy.ones(2)*Value[1]+0.5,
                           lw=0.25,
                           color='black')
            PylabAxis.text(float(XXMax)+XXOffset,
                           float(Value[0]+Value[1])*0.5,
                           r'${\rm '+Key+'}$',
                           verticalalignment='center')
        for Entry in XLeft:
            PylabAxis.plot(scipy.array([Entry,Entry]),
                           scipy.array(PylabAxis.get_ylim()),
                           lw=0.25,
                           color='black')
        for Entry in XRight:
            PylabAxis.plot(scipy.array([Entry,Entry]),
                           scipy.array(PylabAxis.get_ylim()),
                           lw=0.25,
                           color='black')

#        PylabAxis.set_ylim([0,164])
        PylabAxis.spines['right'].set_visible(False)
        PylabAxis.spines['top'].set_visible(False)
        PylabAxis.xaxis.set_ticks_position('bottom')
        PylabAxis.yaxis.set_ticks_position('left')
        PylabAxis.xaxis.set_ticks(XTicks)
        PylabAxis.xaxis.set_ticklabels(XTickLabels)
        for Label in PylabAxis.xaxis.get_ticklabels():
            Label.set_rotation(90)
        pylab.savefig(PlotFile)
    else:
        print '!! NOT IMPLEMENTED YET !!'

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