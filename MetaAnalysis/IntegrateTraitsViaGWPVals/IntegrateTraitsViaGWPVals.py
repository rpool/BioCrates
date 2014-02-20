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

if(False):
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
    MtbNameDict  = {}
    fr = open('Data/Biocrates_Metabolites.csv','r')
    fr.readline()
    for Line in fr:
        LSplit = Line.strip().split(',')
        MtbClassDict[LSplit[2]] = LSplit[4]
        MtbNameDict[LSplit[2]]  = LSplit[1]
    fr.close()

    NodePropertyDict = {}
    for Trait in MtbClassDict.iterkeys():
        NodeId = 'Trait:'+Trait
        NodePropertyDict[NodeId]                       = {}
        NodePropertyDict[NodeId]['NodeType']           = 'Metabolite'
        NodePropertyDict[NodeId]['NodeClass']          = MtbClassDict[Trait]
        NodePropertyDict[NodeId]['NodeName']           = Trait
        NodePropertyDict[NodeId]['NodeBiocratesName']  = MtbNameDict[Trait]
        NodePropertyDict[NodeId]['NodeInENGAGE']       = '1'
        NodePropertyDict[NodeId]['NodeInENGAGEMA']     = '0'

    # Node properties GGM nodes
    for Trait in GGMPartCorr.iterkeys():
        NodeId                                         = 'Trait:'+Trait
        if(not NodePropertyDict.has_key(NodeId)):
            NodePropertyDict[NodeId]                       = {}
            NodePropertyDict[NodeId]['NodeType']           = 'Metabolite'
            NodePropertyDict[NodeId]['NodeClass']          = 'None'
            NodePropertyDict[NodeId]['NodeName']           = 'None'
            NodePropertyDict[NodeId]['NodeBiocratesName']  = 'None'
            NodePropertyDict[NodeId]['NodeInENGAGE']       = '0'
            NodePropertyDict[NodeId]['NodeInENGAGEMA']     = '0'
            if(MtbClassDict.has_key(Trait)):
                NodePropertyDict[NodeId]['NodeClass']          = MtbClassDict[Trait]
                NodePropertyDict[NodeId]['NodeName']           = Trait
                NodePropertyDict[NodeId]['NodeBiocratesName']  = MtbNameDict[Trait]

    # GGM edge properties
    EdgePropertyDict = {}
    for i in xrange(len(Header)-1):
        for j in xrange(i+1,len(Header)):
            EdgeId                   = 'Trait:'+Header[i]+'*'+'Trait:'+Header[j]
            EdgePropertyDict[EdgeId] = {}

            EdgePropertyDict[EdgeId]['GGMEdgeType']  = 'ggm-ggm'
            EdgePropertyDict[EdgeId]['GGMEdgeValue'] = GGMPartCorr[Header[i]][Header[j]]

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
        NodeId = 'Trait:'+Trait
        NodePropertyDict[NodeId]['NodeInENGAGEMA']    = '1'

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
            EdgeId = 'Trait:'+JaccardTraits[i]+'*'+'Trait:'+JaccardTraits[j]
            if(not EdgePropertyDict.has_key(EdgeId)):
                EdgePropertyDict[EdgeId] = {}
            EdgePropertyDict[EdgeId]['JaccardEdgeType']  = 'jaccard-jaccard'
            EdgePropertyDict[EdgeId]['JaccardEdgeValue'] = JaccardDict[JaccardTraits[i]][JaccardTraits[j]]

    # Gene based p-values
    CurrentKeys = NodePropertyDict.keys()
    fr = open('Data/UniqGenes.tsv','r')
    for Line in fr:
        NodeId                                        = 'Gene:'+Line.strip()
        NodePropertyDict[NodeId]                      = {}
        NodePropertyDict[NodeId]['NodeType']          = 'Gene'
        NodePropertyDict[NodeId]['NodeClass']         = 'GeneSymbol'
        NodePropertyDict[NodeId]['NodeName']          = re.sub('Gene:','',NodeId)
        NodePropertyDict[NodeId]['NodeBiocratesName'] = 'None'
        NodePropertyDict[NodeId]['NodeInENGAGE']      = '0'
        NodePropertyDict[NodeId]['NodeInENGAGEMA']    = '0'
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
            NodeId       = 'Trait:'+TraitId
            NodePropertyDict[NodeId]['NodeGroupDeepSplit'+str(s)] = ClusterIndex
            for GId in DataDict[TraitId]['GeneSetAtAlpha_'+str(OptimalAlpha)]:
                NodePropertyDict['Gene:'+GId]['NodeGroupDeepSplit'+str(s)] = ClusterIndex
                NodePropertyDict['Gene:'+GId]['NodeInENGAGEMA']            = '1'
        fr.close()

    fw = open('Data/Nodes.tsv','w')
    fw.write('NodeId\tNodeName\tNodeBiocratesName\tNodeType\tNodeClass\tNodeInENGAGE\tNodeInENGAGEMA\tNodeClusterDeepSplit0\tNodeClusterDeepSplit4\n')
    for Node in NodePropertyDict.iterkeys():
        StringList = []
        StringList.append(Node)
        try:
            StringList.append(NodePropertyDict[Node]['NodeName'])
        except:
            StringList.append('None')
        try:
            StringList.append(NodePropertyDict[Node]['NodeBiocratesName'])
        except:
            StringList.append('None')
        try:
            StringList.append(NodePropertyDict[Node]['NodeType'])
        except:
            StringList.append('None')
        try:
            StringList.append(NodePropertyDict[Node]['NodeClass'])
        except:
            StringList.append('None')
        try:
            StringList.append(NodePropertyDict[Node]['NodeInENGAGE'])
        except:
            StringList.append('0')
        try:
            StringList.append(NodePropertyDict[Node]['NodeInENGAGEMA'])
        except:
            StringList.append('0')
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
            EdgeId = 'Trait:'+Trait+'*'+'Gene:'+G
            pPVal  = -scipy.log10(PVals[scipy.where(Genes==G)[0]][0])
            if(not EdgePropertyDict.has_key(EdgeId)):
                EdgePropertyDict[EdgeId] = {}
            EdgePropertyDict[EdgeId]['TraitGeneEdgeType']  = 'trait-gene'
            EdgePropertyDict[EdgeId]['TraitGeneEdgeValue'] = str(pPVal)

    fw = open('Data/Edges.tsv','w')
    fw.write('SourceId\tTargetId\tSourceName\tTargetName\tEdgeType\tEdgeValue\tGGMEdgeValue\tJaccardEdgeValue\tTraitGeneEdgeValue\n')
    for Edge in EdgePropertyDict.iterkeys():
        StringList = Edge.split('*')
        StringList.extend(re.sub('Gene:','',re.sub('Trait:','',Edge)).split('*'))
        for s0 in ['GGM','Jaccard','TraitGene']:
            try:
                StringList.append(EdgePropertyDict[Edge][s0+'EdgeType'])
                StringList.append(EdgePropertyDict[Edge][s0+'EdgeValue'])
                for s1 in ['GGM','Jaccard','TraitGene']:
                    if(s1==s0):
                        try:
                            StringList.append(EdgePropertyDict[Edge][s1+'EdgeValue'])
                        except:
                            StringList.append('None')
                    else:
                        StringList.append('None')
                fw.write('\t'.join(StringList)+'\n')
                StringList = Edge.split('*')
                StringList.extend(re.sub('Gene:','',re.sub('Trait:','',Edge)).split('*'))
            except:
                StringList = Edge.split('*')
                StringList.extend(re.sub('Gene:','',re.sub('Trait:','',Edge)).split('*'))
    fw.close()