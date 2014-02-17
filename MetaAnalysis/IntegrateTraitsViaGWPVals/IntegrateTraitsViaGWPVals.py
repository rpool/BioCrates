#! /usr/bin/env python

import sys
import os
import scipy
import statsmodels.stats.multitest
import fnmatch
import json
import re

# os.system('lbzip2 -d -k -f /home/r.pool/Work/GWABioCrates/GeneWisePValues/Biocrates_ENGAGE_genewise_pvalues_corrected/PC_aa_C38_4.npy.bz2')
# PVals = scipy.load('/home/r.pool/Work/GWABioCrates/GeneWisePValues/Biocrates_ENGAGE_genewise_pvalues_corrected/PC_aa_C38_4.npy')[1,1:].astype(float)
# Genes = scipy.load('/home/r.pool/Work/GWABioCrates/GeneWisePValues/Biocrates_ENGAGE_genewise_pvalues_corrected/PC_aa_C38_4.npy')[0,1:]
# os.remove('/home/r.pool/Work/GWABioCrates/GeneWisePValues/Biocrates_ENGAGE_genewise_pvalues_corrected/PC_aa_C38_4.npy')

if(False):
    GeneWisePValuePath = '/home/r.pool/Work/GWABioCrates/GeneWisePValues/Biocrates_ENGAGE_genewise_pvalues_corrected'
    GWPValFiles        = os.listdir(GeneWisePValuePath)
    GWPValFiles        = fnmatch.filter(GWPValFiles,'*.npy.bz2')

    DataDict              = {}

    ALvls = scipy.linspace(0.0,0.001,10,False)
    ALvls = scipy.append(ALvls,scipy.linspace(0.001,0.01,9,False))
    ALvls = scipy.append(ALvls,scipy.linspace(0.01,0.05,8,False))
    ALvls = scipy.append(ALvls,scipy.linspace(0.05,0.1,8,False))
    ALvls = scipy.append(ALvls,scipy.linspace(0.1,0.25,16,False))
    ALvls = scipy.append(ALvls,scipy.linspace(0.25,0.75,21,True))
    ALvls = scipy.append(ALvls,scipy.array([0.241]))
    ALvls = scipy.sort(ALvls)

    AlphaLvls = []
    for a in ALvls:
        AlphaLvls.append(a)
    DataDict['AlphaLvls'] = AlphaLvls

    AllGenes = scipy.array([])
    for F in GWPValFiles:
        File            = os.path.join(GeneWisePValuePath,F)
        DecomprFile     = File[:-4]
        Trait           = re.sub('_','.',os.path.basename(DecomprFile)[:-4])
        DataDict[Trait] = {}
        os.system('lbzip2 -d -k -f '+File)
        Data                      = scipy.load(DecomprFile)
        PVals                     = Data[1,1:].astype(float)
        Genes                     = Data[0,1:]
        AllGenes                  = scipy.append(AllGenes,Genes)
        AllGenes                  = scipy.unique(AllGenes)
        NGenes                    = len(Genes)
        DataDict[Trait]['NGenes'] = NGenes
        os.remove(DecomprFile)
        for Alpha in AlphaLvls:
            BH    = statsmodels.stats.multitest.multipletests(pvals=PVals,
                                                              alpha=Alpha,
                                                              method='fdr_bh',
                                                              returnsorted=False)
            DataDict[Trait]['Alpha_'+str(Alpha)]     = Alpha
            DataDict[Trait]['AlphaBonf_'+str(Alpha)] = BH[3]
            BHPVals  = scipy.array(['BHp_value_alpha='+str(Alpha)])
            BHPVals  = scipy.append(BHPVals,BH[1].astype(str))
            Data     = scipy.vstack((Data,BHPVals))
            BHAccept = scipy.array(['BHAccept_alpha='+str(Alpha)])
            BHAccept = scipy.append(BHAccept,BH[0].astype(str))
            DataDict[Trait]['GeneSetAtAlpha_'+str(Alpha)] = scipy.compress(condition=BH[0],
                                                                           a=Genes).tolist()
            Data     = scipy.vstack((Data,BHAccept))
        OutFile  = os.path.join('Data',os.path.basename(DecomprFile))
        scipy.save(file=OutFile,
                   arr=Data)
        os.system('lbzip2 -f '+OutFile)
        print OutFile
    fw = open('Data/DataDict.json','w')
    json.dump(obj=DataDict,
              fp=fw)
    fw.close()
    os.system('lbzip2 -f Data/DataDict.json')

    AllGenesFile = 'Data/UniqGenesOverAllTraits.npy'
    scipy.save(file=AllGenesFile,
               arr=AllGenes)
    os.system('lbzip2 -f '+AllGenesFile)

if(True):
    AllGenesFile = 'Data/UniqGenesOverAllTraits.npy'
    os.system('lbzip2 -d -k '+AllGenesFile+'.bz2')
    AllGenesInGWPValueFiles = scipy.load('Data/UniqGenesOverAllTraits.npy')
    os.remove(AllGenesFile)
#    Output files for CytoScape
    scipy.savetxt(fname='Data/UniqGenes.tsv',
                  X=AllGenesInGWPValueFiles,
                  fmt='%s')
#   Get GO terms of background GeneSet
#     OBORdr = goatools.obo_parser.OBOReader(obo_file='/home/r.pool/Downloads/goatools/gene_ontology.1_2.obo')
#     GODag  = goatools.obo_parser.GODag(obo_file='/home/r.pool/Downloads/goatools/gene_ontology.1_2.obo')
#     print GODag.paths_to_top(term='GO:0003682',
#                              verbose=True)
#     print GODag.viewkeys()
# #     print GODag
# #     for Rec in OBORdr:
# #         print Rec
#
# #   Get OMIM terms of background GeneSet

if(False):
    JsonFile = 'Data/DataDict.json.bz2'
    os.system('lbzip2 -d -f -k '+JsonFile)
    DecomprJsonFile = 'Data/DataDict.json'
    fr              = open(DecomprJsonFile,'r')
    DataDict        = json.load(fp=fr)
    fr.close()
    os.remove(DecomprJsonFile)
    Traits = DataDict.keys()
    del Traits[Traits.index('AlphaLvls')]
    Traits.sort()
    fw = open('Data/IndicesOfTraits.csv','w')
    fw.write('Index,Trait\n')
    for i in xrange(len(Traits)):
        fw.write(str(i)+','+Traits[i]+'\n')
    fw.close()

    MAResultsFile   = 'Data/NewStuffEditAddNearestReportedSNP_PropVar_rMZ_PropH2_16092013_SheetAllData.csv'
    fr              = open(MAResultsFile,'r')
    MAResultsHeader = fr.readline().strip().split(',')
    MATraits        = []
    MAPVals         = []
    for Line in fr:
        LSplit = Line.strip().split(',')
        MATraits.append(LSplit[MAResultsHeader.index('Trait')])
        MAPVals.append(LSplit[MAResultsHeader.index('P-value')])
    fr.close()
    MATraits = scipy.array(MATraits)
    MAPVals  = scipy.array(MAPVals).astype(float)

    GWAlpha   = 5.0e-8
    GWMWAlpha = GWAlpha/46.0

    GWSignTraits   = scipy.unique(MATraits[scipy.where(MAPVals<GWAlpha)[0]])
    GWMWSignTraits = scipy.unique(MATraits[scipy.where(MAPVals<GWMWAlpha)[0]])

    fw = open('Data/NTotalGenesAndGeneSetInfoOfAlpha.tsv','w')
    for Alpha in DataDict['AlphaLvls']:
        if(Alpha==0.0):
            continue
        NTotalGenesOfAlpha   = 0
        UniqGenes            = []
        TraitSetAtAlpha      = []
        for i in xrange(len(Traits)):
            GeneSetAtAlpha      = DataDict[Traits[i]]['GeneSetAtAlpha_'+str(Alpha)]
            NTotalGenesOfAlpha += len(GeneSetAtAlpha)
            UniqGenes.extend(GeneSetAtAlpha)
            if(len(GeneSetAtAlpha)>0):
                TraitSetAtAlpha.append(Traits[i])
        TraitSetAtAlpha  = scipy.array(TraitSetAtAlpha)
        GWIntersection   = scipy.intersect1d(ar1=TraitSetAtAlpha,
                                             ar2=GWSignTraits,
                                             assume_unique=False)
        GWMWIntersection = scipy.intersect1d(ar1=TraitSetAtAlpha,
                                             ar2=GWMWSignTraits,
                                             assume_unique=False)
        GWUnion          = scipy.union1d(ar1=TraitSetAtAlpha,
                                         ar2=GWSignTraits)
        GWMWUnion        = scipy.union1d(ar1=TraitSetAtAlpha,
                                         ar2=GWMWSignTraits)
        fw.write(str(Alpha)+'\t'+\
                 str(NTotalGenesOfAlpha)+'\t'+\
                 str(len(scipy.unique(scipy.array(UniqGenes))))+'\t'+\
                 str(len(TraitSetAtAlpha))+'\t'+\
                 str(len(GWSignTraits))+'\t'+\
                 str(len(GWIntersection))+'\t'+\
                 str(float(len(GWIntersection))/float(len(GWUnion)))+'\t'+\
                 str(len(TraitSetAtAlpha))+'\t'+\
                 str(len(GWMWSignTraits))+'\t'+\
                 str(len(GWMWIntersection))+'\t'+\
                 str(float(len(GWMWIntersection))/float(len(GWMWUnion)))+'\n')
    fw.close()

if(False):
    JsonFile = 'Data/DataDict.json.bz2'
    os.system('lbzip2 -d -f -k '+JsonFile)
    DecomprJsonFile = 'Data/DataDict.json'
    fr              = open(DecomprJsonFile,'r')
    DataDict        = json.load(fp=fr)
    fr.close()
    os.remove(DecomprJsonFile)
    Traits = DataDict.keys()
    del Traits[Traits.index('AlphaLvls')]
    Traits.sort()
    fw = open('Data/IndicesOfTraits.csv','w')
    fw.write('Index,Trait\n')
    for i in xrange(len(Traits)):
        fw.write(str(i)+','+Traits[i]+'\n')
    fw.close()

    OptimalAlpha = 0.241
#     for Alpha in DataDict['AlphaLvls']:
    for Alpha in [OptimalAlpha]:
        JaccardArray = scipy.array([[0.0]*len(Traits)]*len(Traits))
        fw = open('Data/JaccardArrayAlpha'+str(Alpha)+'.csv','w')
        Indices = []
        for i in xrange(len(Traits)):
#             print Traits[i],len(DataDict[Traits[i]]['GeneSetAtAlpha_'+str(Alpha)])
            if(len(DataDict[Traits[i]]['GeneSetAtAlpha_'+str(Alpha)])==0):
                continue
            Indices.append(i)
        fw.write(','.join(scipy.array(Traits)[Indices].tolist())+'\n')
        for i in Indices:
            StringList = []
            for j in Indices:
                GSet_i            = scipy.array(DataDict[Traits[i]]['GeneSetAtAlpha_'+str(Alpha)])
                GSet_j            = scipy.array(DataDict[Traits[j]]['GeneSetAtAlpha_'+str(Alpha)])
                Jaccard           = float(len(scipy.intersect1d(GSet_i,GSet_j)))
                Jaccard          /= max(1.0e-10,float(len(scipy.union1d(GSet_i,GSet_j))))
                JaccardArray[i,j] = Jaccard
                StringList.append(str(JaccardArray[i,j]))
            JaccardArray[i,i] = 1.0
            fw.write(','.join(StringList)+'\n')
        fw.close()

        # Yoshiko clustering
        Thresholds = [0.0,0.01,0.02,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4]
        for t in Thresholds:
            YoshikoThresholdFile = 'Data/Alpha_'+str(Alpha)+'_Threshold_'+str(t)+'.yok'
            fw = open(YoshikoThresholdFile,'w')
            fw.write(str(len(Traits))+'\n')
            for T in Traits:
                fw.write(T+'\n')
            for i in xrange(len(Traits)-1):
                fw.write(' '.join((JaccardArray[i,i+1:]-t).astype(str).tolist())+'\n')
            fw.close()
            YoshikoOutFile = re.sub('.yok','.o',YoshikoThresholdFile)
            YoshikoErrFile = re.sub('.yok','.e',YoshikoThresholdFile)
            YoshikoSolFile = re.sub('.yok','.sol',YoshikoThresholdFile)
            os.system('yoshiko2.0 -F 0 -f '+YoshikoThresholdFile+' -o '+YoshikoSolFile+' -v 5 -r 11111 -m 10 > '+YoshikoOutFile+' 2> '+YoshikoErrFile)
            fw           = open('Data/ClusterSummary_Alpha_'+str(Alpha)+'_Threshold_'+str(t)+'.tsv','w')
            HeaderString = 'ClusterIndex\tTraits\tNTraits\tGeneSetUnion\tNGenesInUnion'
            fw.write(HeaderString+'\n')
            fr        = open(YoshikoSolFile+'_0.csv','r')
            NClusters = 0
            for Line in fr:
                NClusters += 1
                LSplit     = Line.strip().split()
                LSplit.sort()
                fw.write(str(NClusters)+'\t')
                fw.write(str(LSplit)+'\t')
                fw.write(str(len(LSplit))+'\t')
                GeneSetUnion = []
                for Trait in LSplit:
                    GeneSetUnion.extend(DataDict[Trait]['GeneSetAtAlpha_'+str(Alpha)])
                GeneSetUnion.sort()
                fw.write(str(GeneSetUnion)+'\t')
                fw.write(str(len(GeneSetUnion)))
                fw.write('\n')
            fr.close()
            fw.close()

        # R dynamic tree cut clustering
        os.system('R CMD BATCH ./CutTree.R') # yields two files for different values of 'deepSplit'
        for s in [0,4]:
            DeepSplitFile = 'Data/TreeCutJaccardArrayAlpha'+str(OptimalAlpha)+'DeepSplit'+str(s)+'.csv'

            fr = open(DeepSplitFile,'r')
            fw = open('Data/TreeCutClusterSummaryAlpha'+str(OptimalAlpha)+'DeepSplit'+str(s)+'.tsv','w')
            fr.readline() # skip header
            Clusters  = []
            ClustDict = {}
            for Line in fr:
                LSplit       = Line.strip().split(',')
                ClusterIndex = int(LSplit[-1])
                Clusters.append(ClusterIndex)
                if(not ClustDict.has_key(ClusterIndex)):
                    ClustDict[ClusterIndex]           = {}
                    ClustDict[ClusterIndex]['Traits'] = []
                    ClustDict[ClusterIndex]['Genes']  = []
                ClustDict[ClusterIndex]['Traits'].append(re.sub('\"','',LSplit[0]))
                ClustDict[ClusterIndex]['Genes'].extend(DataDict[re.sub('\"','',LSplit[0])]['GeneSetAtAlpha_'+str(OptimalAlpha)])
            Clusters = list(set(Clusters))
            Clusters.sort()
            fw.write('ClusterIndex\tTraits\tNTraits\tGenes\tNGenes\n')
            for c in Clusters:
                fw.write(str(c)+'\t'+\
                         '|'.join(ClustDict[c]['Traits'])+'\t'+\
                         str(len(ClustDict[c]['Traits']))+'\t'+\
                         '|'.join(list(set(ClustDict[c]['Genes'])))+'\t'+\
                         str(len(list(set(ClustDict[c]['Genes']))))+'\n')
            fr.close()
            fw.close()

if(True):
    # Read GGM
    fr     = open('Data/ggm_cn_sheet_partial.csv','r')
    Header = fr.readline().strip().split(',')
    del Header[Header.index('')]
    for i in xrange(len(Header)):
        KORATrait       = Header[i]
        ConventionTrait = re.sub('[\(,\),\-,:,\ ,\/]','.',KORATrait)
        while(re.search('\.\.',ConventionTrait)):
            ConventionTrait = re.sub('\.\.','.',ConventionTrait)
        if(re.match('C[0-9].0',ConventionTrait)):
            ConventionTrait = re.sub('\.0','',ConventionTrait)
        if(re.match('C1[0-9].0',ConventionTrait)):
            ConventionTrait = re.sub('\.0','',ConventionTrait)
        Header[i] = ConventionTrait
    GGMPartCorr = {}
    for Trait in Header:
        GGMPartCorr[Trait] = {}
    for Line in fr:
        LSplit = Line.strip().split(',')
        Trait  = LSplit[0]
        Trait  = re.sub('[\(,\),\-,:,\ ,\/]','.',Trait)
        while(re.search('\.\.',Trait)):
            Trait = re.sub('\.\.','.',Trait)
        if(re.match('C[0-9].0',Trait)):
            Trait = re.sub('\.0','',Trait)
        if(re.match('C1[0-9].0',Trait)):
            Trait = re.sub('\.0','',Trait)
        for i in xrange(len(Header)):
            GGMPartCorr[Header[i]][Trait] = LSplit[i+1]
    fr.close()

    # Metabolite classes
    MtbClassDict = {}
    fr = open('Data/Biocrates_Metabolites.csv','r')
    fr.readline()
    for Line in fr:
        LSplit = Line.strip().split(',')
        MtbClassDict[LSplit[2]] = LSplit[4]
    fr.close()

    # Node properties GGM nodes
    NodePropertyDict = {}
    for Trait in GGMPartCorr.iterkeys():
        NodePropertyDict[Trait]              = {}
        NodePropertyDict[Trait]['NodeType']  = 'Metabolite'
        NodePropertyDict[Trait]['NodeClass'] = 'None'
        if(MtbClassDict.has_key(Trait)):
            NodePropertyDict[Trait]['NodeClass'] = MtbClassDict[Trait]

    # GGM edge properties
    EdgePropertyDict = {}
    for i in xrange(len(Header)-1):
        for j in xrange(i+1,len(Header)):
            EdgeId                   = Header[i]+'*'+Header[j]
            EdgePropertyDict[EdgeId] = {}

            EdgePropertyDict[EdgeId]['EdgeType']  = 'ggm-ggm'
            EdgePropertyDict[EdgeId]['EdgeValue'] = GGMPartCorr[Header[i]][Header[j]]

    # Read ENGAGE data
    JsonFile = 'Data/DataDict.json.bz2'
    os.system('lbzip2 -d -f -k '+JsonFile)
    DecomprJsonFile = 'Data/DataDict.json'
    fr              = open(DecomprJsonFile,'r')
    DataDict        = json.load(fp=fr)
    fr.close()
    os.remove(DecomprJsonFile)
    Traits = DataDict.keys()
    del Traits[Traits.index('AlphaLvls')]
    Traits.sort()
    for i in xrange(len(Traits)):
        Traits[i] = str(Traits[i])

    # ENGAGE node properties
    for Trait in Traits:
        if(NodePropertyDict.has_key(Trait)):
            continue
        NodePropertyDict[Trait]              = {}
        NodePropertyDict[Trait]['NodeType']  = 'Metabolite'
        NodePropertyDict[Trait]['NodeClass'] = 'None'
        if(MtbClassDict.has_key(Trait)):
            NodePropertyDict[Trait]['NodeClass'] = MtbClassDict[Trait]

    # ENGAGE edge properties
    OptimalAlpha = 0.241
    fr = open('Data/JaccardArrayAlpha'+str(OptimalAlpha)+'.csv','r')
    JaccardTraits = fr.readline().strip().split(',')
    JaccardDict   = {}
    for Trait in JaccardTraits:
        JaccardDict[str(Trait)] = {}
    Index = 0
    for Line in fr:
        LSplit = Line.strip().split(',')
        for i in xrange(len(JaccardTraits)):
            JaccardDict[str(JaccardTraits[Index])][str(JaccardTraits[i])] = LSplit[i]
        Index += 1
    fr.close()
    for i in xrange(len(JaccardTraits)-1):
        for j in xrange(i+1,len(JaccardTraits)):
            EdgeId = JaccardTraits[i]+'*'+JaccardTraits[j]
            if(not EdgePropertyDict.has_key(EdgeId)):
                EdgePropertyDict[EdgeId] = {}
            EdgePropertyDict[EdgeId]['EdgeType']  = 'jaccard-jaccard'
            EdgePropertyDict[EdgeId]['EdgeValue'] = JaccardDict[JaccardTraits[i]][JaccardTraits[j]]

    # Gene based p-values
    fr = open('Data/UniqGenes.tsv','r')
    for Line in fr:
        NodePropertyDict[Line.strip()]              = {}
        NodePropertyDict[Line.strip()]['NodeType']  = 'Gene'
        NodePropertyDict[Line.strip()]['NodeClass'] = 'GeneSymbol'
    fr.close()

    for s in [0,4]:
        for Node in NodePropertyDict.iterkeys():
            NodePropertyDict[Node]['NodeGroupDeepSplit'+str(s)] = 'None'

        DeepSplitFile = 'Data/TreeCutJaccardArrayAlpha'+str(OptimalAlpha)+'DeepSplit'+str(s)+'.csv'
        fr            = open(DeepSplitFile,'r')
        fr.readline() # skip header
        for Line in fr:
            LSplit       = Line.strip().split(',')
            ClusterIndex = LSplit[-1]
            TraitId      = re.sub('\"','',LSplit[0])
            NodePropertyDict[TraitId]['NodeGroupDeepSplit'+str(s)] = ClusterIndex
            for GId in DataDict[TraitId]['GeneSetAtAlpha_'+str(OptimalAlpha)]:
                NodePropertyDict[GId]['NodeGroupDeepSplit'+str(s)] = ClusterIndex
        fr.close()

    fw = open('Data/Nodes.tsv','w')
    fw.write('NodeId\tNodeType\tNodeClass\tNodeClusterDeepSplit0\tNodeClusterDeepSplit4\n')
    for Node in NodePropertyDict.iterkeys():
        StringList = []
        StringList.append(Node)
        try:
            StringList.append(NodePropertyDict[Node]['NodeType'])
        except:
            StringList.append('None')
        try:
            StringList.append(NodePropertyDict[Node]['NodeClass'])
        except:
            StringList.append('None')
        for s in [0,4]:
            StringList.append(NodePropertyDict[Node]['NodeGroupDeepSplit'+str(s)])
        fw.write('\t'.join(StringList)+'\n')
    fw.close()

    for Trait in JaccardTraits:
        print Trait
        F    = re.sub('\.','_',Trait)+'.npy.bz2'
        File = os.path.join('Data',F)
        os.system('lbzip2 -d -k -f '+File)
        DecomprFile = File[:-4]
        Data        = scipy.load(DecomprFile)
        PVals       = Data[1,1:].astype(float)
        Genes       = Data[0,1:]
        os.remove(DecomprFile)
        GeneSetAtAlpha = DataDict[Trait]['GeneSetAtAlpha_'+str(OptimalAlpha)]
        for G in GeneSetAtAlpha:
            EdgeId = Trait+'*'+G
            pPVal  = -scipy.log10(PVals[scipy.where(Genes==G)[0]][0])
            if(not EdgePropertyDict.has_key(EdgeId)):
                EdgePropertyDict[EdgeId] = {}
            EdgePropertyDict[EdgeId]['EdgeType']  = 'trait-gene'
            EdgePropertyDict[EdgeId]['EdgeValue'] = str(pPVal)

    fw = open('Data/Edges.tsv','w')
    fw.write('Source\tTarget\tEdgeType\tEdgeValue\n')
    for Edge in EdgePropertyDict.iterkeys():
        StringList = Edge.split('*')
        try:
            StringList.append(EdgePropertyDict[Edge]['EdgeType'])
        except:
            StringList.append('None')
        try:
            StringList.append(EdgePropertyDict[Edge]['EdgeValue'])
        except:
            StringList.append('None')
        fw.write('\t'.join(StringList)+'\n')
    fw.close()

#         Clusters = scipy.cluster.hierarchy.fclusterdata(X=JaccardArray,
#                                                         t=0.8,
#                                                         criterion='inconsistent')
#         print len(Clusters),len(scipy.unique(Clusters))
#         print Clusters

# if(True):
#     # Filter PVals on overlap w/ PINA
#     PINAGenes = scipy.genfromtxt(fname='/home/r.pool/workspace/BioCrates/GGM/CWIGGMCytoscape/PINA/PINAOnlyHGNCNoSUMONoUBCNonRedundant.tsv',
#                                  dtype=str,
#                                  unpack=True,
#                                  skip_header=1)
#     PINAG     = networkx.Graph()
#     for i in range(len(PINAGenes[0])):
#         PINAG.add_edge(PINAGenes[0,i],PINAGenes[1,i])
#     PINAG.remove_edges_from(PINAG.selfloop_edges())
#     PINAGenes  = scipy.unique(scipy.append(PINAGenes[0],PINAGenes[1])).tolist()
#     GList      = networkx.connected_component_subgraphs(PINAG)
#
#     IndexOfLargestG = 0
#     SizeOfLargestG  = 0
#     for i in xrange(len(GList)):
#         if(GList[i].size()>SizeOfLargestG):
#             SizeOfLargestG  = GList[i].size()
#             IndexOfLargestG = i
#     PINAG = GList[IndexOfLargestG]
#
#     KeepFilter = [False]*len(Genes)
#     for i in xrange(len(Genes)):
#         if(Genes[i] in PINAGenes):
#             KeepFilter[i] = True
#         else:
#             try:
#                 PINAG.remove_node(Genes[i])
#             except:
#                 pass
#     KeepFilter = scipy.array(KeepFilter)
#     PVals = scipy.compress(KeepFilter,PVals)
#     Genes = scipy.compress(KeepFilter,Genes)
#
# #Indices = scipy.where(PVals<1.0e-1)[0]
# #PVals   = PVals[Indices]
# #Genes   = Genes[Indices]
#
# BH = statsmodels.stats.multitest.multipletests(pvals=PVals,
#                                                alpha=0.05,
#                                                method='fdr_bh',
#                                                returnsorted=False)
# #print BH
#
# scipy.savetxt(fname='PythonBHCorrected.txt',
#               X=BH[1],
#               fmt='%10.10e')
#
# ChiSq = scipy.stats.chi2.isf(PVals,1)
#
# scipy.savetxt(fname='Chi2.txt',
#               X=ChiSq,
#               fmt='%10.10e')
#
# PValsRvs = scipy.stats.chi2.sf(scipy.stats.chi2.rvs(1,size=len(PVals)),1)
# scipy.savetxt(fname='CheckDistribution.txt',
#               X=scipy.array([scipy.sort(-scipy.log10(PValsRvs)),scipy.sort(-scipy.log10(PVals))]).T,
#               fmt='%10.10e')
# #sys.exit()
#
#
# FChiSq = scipy.stats.chi2.pdf(ChiSq,1)
#
# scipy.savetxt(fname='FChi2.txt',
#               X=scipy.array([ChiSq,FChiSq,PVals]).T,
#               fmt='%10.10e')
#
# fw = open('PValsCheck.txt','w')
# for i in xrange(len(Genes)):
#     fw.write(Genes[i]+' '+str(PVals[i])+'\n')
# fw.close()
#
# Alpha        = 0.05
# NIndTests    = 46
# NGenes       = len(PVals)
# #BonfSignLvl  = Alpha/float(NGenes*NIndTests)
# #BonfSignLvl  = Alpha/float(NGenes)
# BonfSignLvl  = Alpha/float(NGenes)
# FBonfSignLvl = scipy.stats.chi2.pdf(scipy.stats.chi2.isf(BonfSignLvl,1),1)
# #print BonfSignLvl, FBonfSignLvl
# #print BonfSignL)l*float(NIndTests), FBonfSignLvl
# #print len(scipy.where(PVals<(BonfSignLvl*float(NIndTests)))[0])
#
# scipy.savetxt(fname='FChi2Bonferroni.txt',
#               X=scipy.array([ChiSq,FChiSq,[FBonfSignLvl]*NGenes,PVals,[BonfSignLvl]*NGenes]).T,
#               fmt='%10.10e')
#
# Header   = ['GeneSymbolHGNC',
#             'PVal',
#             'BonferroniAlpha0.01',
#             'BonferroniAlpha0.05',
#             'BonferroniAlpha0.1',
#             'WeightBonferroniAlpha0.01',
#             'WeightBonferroniAlpha0.05',
#             'WeightBonferroniAlpha0.1',
#             'ChiSq',
#             'PdfChiSq',
#             'BonferroniChiSqCut.01',
#             'BonferroniChiSqCut0.05',
#             'BonferroniChiSqCut0.1',
#             'WeightBonferroniChiSqCut.01',
#             'WeightBonferroniChiSqCut0.05',
#             'WeightBonferroniChiSqCut0.1']
# scipy.set_printoptions(formatter={'float': lambda x: format(x,'12.10e')})
# Formats  = ['S100']
# Formats.extend(['f8']*(len(Header)-1))
# OutArray = scipy.empty(shape=NGenes,dtype={'names':Header,'formats':Formats})
# OutArray['GeneSymbolHGNC']              = Genes
# OutArray['PVal']                        = PVals
# OutArray['BonferroniAlpha0.01']         = scipy.array([0.01]*NGenes)/float(NGenes)
# OutArray['BonferroniAlpha0.05']         = scipy.array([0.05]*NGenes)/float(NGenes)
# OutArray['BonferroniAlpha0.1']          = scipy.array([0.1]*NGenes)/float(NGenes)
# OutArray['WeightBonferroniAlpha0.01']   = -scipy.log10(PVals)+scipy.log10(scipy.array([0.01]*NGenes)/float(NGenes))
# OutArray['WeightBonferroniAlpha0.05']   = -scipy.log10(PVals)+scipy.log10(scipy.array([0.05]*NGenes)/float(NGenes))
# OutArray['WeightBonferroniAlpha0.1']    = -scipy.log10(PVals)+scipy.log10(scipy.array([0.1]*NGenes)/float(NGenes))
# OutArray['ChiSq']                       = ChiSq
# OutArray['PdfChiSq']                    = FChiSq
# OutArray['BonferroniChiSqCut.01']       = scipy.stats.chi2.pdf(scipy.stats.chi2.isf(scipy.array([0.01]*NGenes)/float(NGenes),1),1)
# OutArray['BonferroniChiSqCut0.05']      = scipy.stats.chi2.pdf(scipy.stats.chi2.isf(scipy.array([0.05]*NGenes)/float(NGenes),1),1)
# OutArray['BonferroniChiSqCut0.1']       = scipy.stats.chi2.pdf(scipy.stats.chi2.isf(scipy.array([0.1]*NGenes)/float(NGenes),1),1)
# OutArray['WeightBonferroniChiSqCut.01'] = -scipy.log10(FChiSq)+scipy.log10(scipy.stats.chi2.pdf(scipy.stats.chi2.isf(scipy.array([0.01]*NGenes)/float(NGenes),1),1))
# OutArray['WeightBonferroniChiSqCut0.05']= -scipy.log10(FChiSq)+scipy.log10(scipy.stats.chi2.pdf(scipy.stats.chi2.isf(scipy.array([0.05]*NGenes)/float(NGenes),1),1))
# OutArray['WeightBonferroniChiSqCut0.1'] = -scipy.log10(FChiSq)+scipy.log10(scipy.stats.chi2.pdf(scipy.stats.chi2.isf(scipy.array([0.1]*NGenes)/float(NGenes),1),1))
#
# CommentString = '# This file contains the p-values, pdf_chi2(p-values) and their associated Bonferroni-corrected weights (NGenes='+str(NGenes)+')'+'\n'+\
#                 '# at Alpa levels [0.01,0.05,0.1] for metabolite \"PC_aa_C38_4\".'+'\n'+\
#                 '# \n'+\
#                 '# Column names are on the next row.\n'+\
#                 '\t'.join(Header)+'\n'
# Format = ['%s']
# Format.extend(['%12.10e']*(len(Header)-1))
# fw = open('Weights_PC_aa_C38_4.tsv','w')
# fw.write(CommentString)
# scipy.savetxt(fname=fw,
#               X=OutArray,
#               fmt=Format,
#               delimiter='\t')
# fw.close()
#
#
# #scipy.savetxt(fname='Weights_PC_aa_C38_4.tsv',
# #              X=OutArray,
# #              fmt='%s',
# #              delimiter='\t',
# ##              delimiter='\t')
# ##              header='\t'.join(Header))
# #              header=CommentString+'\n'+'\t'.join(Header))
#
# #sys.exit()
#
# print len(scipy.where(-scipy.log10(PVals/BonfSignLvl)>=0.0)[0])
# print len(scipy.where(-scipy.log10(FChiSq/FBonfSignLvl)>=0.0)[0])
# print len(scipy.where(BH[0])[0])
# print Genes[scipy.where(-scipy.log10(PVals/BonfSignLvl)>=0.0)[0]]
# print Genes[scipy.where(-scipy.log10(FChiSq/FBonfSignLvl)>=0.0)[0]]
# print Genes[scipy.where(BH[0])[0]]
# print PVals[scipy.where(BH[0])[0]]
#
# SelectedGenes        = Genes[scipy.where(-scipy.log10(PVals/BonfSignLvl)>=0.0)[0]]
# SelectionG           = networkx.Graph()
# for i in xrange(len(SelectedGenes)-1):
#     for j in xrange(i+1,len(SelectedGenes)):
#         print SelectedGenes[i],SelectedGenes[j]
#         SelectionG.add_path(networkx.dijkstra_path(PINAG,SelectedGenes[i],SelectedGenes[j],'distance'))
#         print '**',networkx.dijkstra_path(PINAG,SelectedGenes[i],SelectedGenes[j],'distance')
# SelectionG.remove_edges_from(SelectionG.selfloop_edges())
# SelectionG = networkx.Graph(SelectionG)
#
# for i in xrange(1):
#     for N in SelectionG.nodes():
#         for NN in PINAG.neighbors(N):
#             if(NN in Genes.tolist()):
#                 SelectionG.add_edge(N,NN)
#     SelectionG.remove_edges_from(SelectionG.selfloop_edges())
#     SelectionG = networkx.Graph(SelectionG)
# #print SelectionG.size()
# #print PINAG.size()
#
# ArgSort     = scipy.argsort(PVals)
# SortedPVals = scipy.sort(PVals)
# RankArray   = scipy.arange(1,len(PVals)+1,dtype=float)
# ReverseRankArray = RankArray[::-1]
# #print
# #print len(PVals)
# #print SortedPVals
# #print RankArray
# #print ReverseRankArray
# #print RankArray[-1]/ReverseRankArray
# #print SortedPVals*(RankArray[-1]/RankArray)
# #print scipy.sort(BH[1])
# scipy.savetxt(fname='PythonBHCorrectedHomeBrew.txt',
#               X=SortedPVals*(RankArray[-1]/RankArray),
#               fmt='%10.10e')
# scipy.savetxt(fname='PythonBHCorrectedHomeBrewViaLn.txt',
#               X=scipy.exp(scipy.log(SortedPVals)+scipy.log(RankArray[-1])-scipy.log(RankArray)),
#               fmt='%10.10e')
# scipy.savetxt(fname='PythonBHCorrectedHomeBrewLn.txt',
#               X=-(scipy.log10(SortedPVals)+scipy.log10(RankArray[-1])-scipy.log10(RankArray)),
#               fmt='%10.10e')
#
# #print len(scipy.where(SortedPVals*(RankArray[-1]/RankArray)<=Alpha)[0])
# SortedGenes     = Genes[ArgSort]
# SelectedIndices = []
# for N in SelectionG.nodes():
# #    print N,
#     SelectedIndices.append(scipy.where(SortedGenes==N)[0][0])
# #    print SelectedIndices[-1]
#
# SelectedIndices = scipy.unique(scipy.array(SelectedIndices))
# #print SortedGenes[scipy.where((SortedPVals*(RankArray[-1]/RankArray))<=Alpha)[0]]
# #print SortedPVals[scipy.where((SortedPVals*(RankArray[-1]/RankArray))<=Alpha)[0]]
# #print SortedGenes[scipy.where((SortedPVals*(RankArray[-1]/RankArray))<=(Alpha*2.0))[0]]
# #print SortedPVals[scipy.where((SortedPVals*(RankArray[-1]/RankArray))<=(Alpha*2.0))[0]]
#
# scipy.savetxt(fname='NeighborSelection.txt',
#               X=scipy.array([scipy.arange(0,len(SelectedIndices),dtype=float),
#                              SortedPVals[SelectedIndices],
#                              scipy.array(SortedPVals*(RankArray[-1]/RankArray))[SelectedIndices]]).T,
#               fmt='%10.10e')
