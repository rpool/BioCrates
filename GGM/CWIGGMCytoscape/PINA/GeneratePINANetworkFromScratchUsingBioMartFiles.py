#! /usr/bin/env python
import sys
import os
import scipy
import numpy
import itertools

# Generate the PINA network from scratch:

# Read in the PINA data (file "PINA_HSA-20121210.sif.txt", note that these are Uniprot IDs or ACs):
PINAUniprot = scipy.genfromtxt(fname='PINA_HSA-20121210.sif.txt',
                               dtype=str,
                               delimiter=' ',
                               usecols=(0,2),
                               unpack=True)

# Read in the mapping table downloaded from Ensembl BioMart by querying the all unique UniprotKB Accessions (file "BioMartUniprotAC_or_ID_to_HGNCSymbol.tsv")
fr                            = open('BioMartUniprotAC_or_ID_to_HGNCSymbol.tsv','r')
BioMartUniprot2HGNCSymbolsHdr = fr.readline().strip().split('\t')
fr.close()
BioMartUniprot2HGNCSymbols = scipy.genfromtxt(fname='BioMartUniprotAC_or_ID_to_HGNCSymbol.tsv',
                                              dtype=str,
                                              delimiter='\t',
                                              skip_header=1,
                                              unpack=True)

# Check if all PINA UniprotKB IDs are reported in the BioMart file:
AllPINAUniprotKBIDs          = scipy.unique(PINAUniprot[0])
AllPINAUniprotKBIDs          = scipy.append(AllPINAUniprotKBIDs,PINAUniprot[1])
AllPINAUniprotKBIDs          = scipy.unique(AllPINAUniprotKBIDs)
AllUNIProtKBIDsInBioMart     = scipy.unique(BioMartUniprot2HGNCSymbols[BioMartUniprot2HGNCSymbolsHdr.index('UniProt/SwissProt Accession')])
BioMartUniprotKBIDsNotInPINA = scipy.setdiff1d(ar1=AllPINAUniprotKBIDs,
                                               ar2=AllUNIProtKBIDsInBioMart,
                                               assume_unique=False)
fw = open('BioMartUniprotKBIDsNotInPINA.txt','w')
for i in xrange(len(BioMartUniprotKBIDsNotInPINA)):
    fw.write(BioMartUniprotKBIDsNotInPINA[i]+'\n')
fw.close()
PINAUniprotKBIDsNotInBioMart = scipy.setdiff1d(ar1=AllUNIProtKBIDsInBioMart,
                                               ar2=AllPINAUniprotKBIDs,
                                               assume_unique=False)
fw = open('PINAUniprotKBIDsNotInBioMart.txt','w')
for i in xrange(len(PINAUniprotKBIDsNotInBioMart)):
    fw.write(PINAUniprotKBIDsNotInBioMart[i]+'\n')
fw.close()

sys.exit()

# Map PINA network to HGNC Gene symbols
UniprotAcsInPINA        = scipy.unique(BioMartUniprot2HGNCSymbols[BioMartUniprot2HGNCSymbolsHdr.index('UniProt/SwissProt Accession')])
UniprotAcID2HGNCSymbol  = {}
MultipleHGNCsForUniprot = {}
for i in xrange(len(UniprotAcsInPINA)):
    UniprotAc        = UniprotAcsInPINA[i]
    Indices          = scipy.where(BioMartUniprot2HGNCSymbols[BioMartUniprot2HGNCSymbolsHdr.index('UniProt/SwissProt Accession')]==UniprotAc)[0]
    HGNCSymbol       = scipy.unique(BioMartUniprot2HGNCSymbols[BioMartUniprot2HGNCSymbolsHdr.index('HGNC symbol'),Indices])

    UniprotAcID2HGNCSymbol[UniprotAc] = []
    L = 0
    for H in HGNCSymbol:
        if(H==''):
            continue
        UniprotAcID2HGNCSymbol[UniprotAc].append(H)
        L += 1
    if(L>1):
        MultipleHGNCsForUniprot[UniprotAc] = HGNCSymbol

# Convert UniprotKB Accessions in PINA to HGNC symbols
PINAHGNC = []
fw       = open('PINAHGNC.tsv','w')
fw.write('NodeA_HGNC\tNodeB_HGNC\tNodeAUniprotKB\tNodeBUniprotKB\n')
for i in xrange(max(scipy.shape(PINAUniprot))):
    NodeAUniProt = PINAUniprot[0,i]
    NodeAList    = ['None']
    if(UniprotAcID2HGNCSymbol.has_key(NodeAUniProt)):
        NodeAList = UniprotAcID2HGNCSymbol[NodeAUniProt]
        while('' in NodeAList):
            del NodeAList[NodeAList.index('')]
    NodeBUniProt = PINAUniprot[1,i]
    NodeBList    = ['None']
    if(UniprotAcID2HGNCSymbol.has_key(NodeBUniProt)):
        NodeBList = UniprotAcID2HGNCSymbol[NodeBUniProt]
        while('' in NodeBList):
            del NodeBList[NodeBList.index('')]
    Product = list(itertools.product(NodeAList,NodeBList))
    for j in xrange(len(Product)):
        fw.write(Product[j][0]+'\t'+
                 Product[j][1]+'\t'+
                 NodeAUniProt+'\t'+
                 NodeBUniProt+'\t\n')
fw.close()

# Find unmatched UniprotKBIDs
PINAHGNC = scipy.genfromtxt(fname='PINAHGNC.tsv',
                            dtype=str,
                            delimiter='\t',
                            skip_header=1,
                            unpack=True)
UnmatchedUniprotKBIDs = scipy.array([])
UnmatchedUniprotKBIDs = scipy.append(UnmatchedUniprotKBIDs,PINAHGNC[2,scipy.where(PINAHGNC[0]=='None')[0]])
UnmatchedUniprotKBIDs = scipy.append(UnmatchedUniprotKBIDs,PINAHGNC[3,scipy.where(PINAHGNC[1]=='None')[0]])
UnmatchedUniprotKBIDs = scipy.unique(UnmatchedUniprotKBIDs)
fw = open('UnmatchedUniprotKBIDs.txt','w')
for i in xrange(len(UnmatchedUniprotKBIDs)):
    fw.write(UnmatchedUniprotKBIDs[i]+'\n')
fw.close()

# Read in union of all genes of all gene-wise p-values over all metabolites (files "UniqGeneSymbols.dat", note that the entries are HGNC Gene symbols):
UniqHGNCSymbolsInENGAGEData = scipy.genfromtxt(fname='UniqGeneSymbols.dat',
                                               dtype=str,
                                               delimiter='\t',
                                               skip_header=1,
                                               unpack=True)

# Determine overlap between PINA and ENGAGE set generated by VEGAS:
GeneSymbolsInPINA = scipy.array([])
GeneSymbolsInPINA = scipy.append(GeneSymbolsInPINA,PINAHGNC[0])
GeneSymbolsInPINA = scipy.append(GeneSymbolsInPINA,PINAHGNC[1])
GeneSymbolsInPINA = GeneSymbolsInPINA[scipy.where(GeneSymbolsInPINA!='None')[0]]
GeneSymbolsInPINA = scipy.unique(GeneSymbolsInPINA)

ENGAGEGeneSymbolsNotInPINA = scipy.setdiff1d(ar1=UniqHGNCSymbolsInENGAGEData,
                                             ar2=GeneSymbolsInPINA,
                                             assume_unique=True)
fw = open('ENGAGEGeneSymbolsNotInPINA.txt','w')
for i in xrange(len(ENGAGEGeneSymbolsNotInPINA)):
    fw.write(ENGAGEGeneSymbolsNotInPINA[i]+'\n')
fw.close()

PINAGeneSymbolsNotInENGAGE = scipy.setdiff1d(ar1=GeneSymbolsInPINA,
                                             ar2=UniqHGNCSymbolsInENGAGEData,
                                             assume_unique=True)
fw = open('PINAGeneSymbolsNotInENGAGE.txt','w')
for i in xrange(len(PINAGeneSymbolsNotInENGAGE)):
    fw.write(PINAGeneSymbolsNotInENGAGE[i]+'\n')
fw.close()

# Remove unmatched UniprotKBIDs:
FilterArray  = (PINAHGNC[0]!='None')
FilterArray *= (PINAHGNC[1]!='None')

# Remove Ubiquitin-related parts from the complete PINA network:
FilterArray *= (PINAHGNC[0]!='UBC')
FilterArray *= (PINAHGNC[0]!='SUMO1')
FilterArray *= (PINAHGNC[0]!='SUMO2')
FilterArray *= (PINAHGNC[0]!='SUMO3')
FilterArray *= (PINAHGNC[0]!='SUMO4')
FilterArray *= (PINAHGNC[1]!='UBC')
FilterArray *= (PINAHGNC[1]!='SUMO1')
FilterArray *= (PINAHGNC[1]!='SUMO2')
FilterArray *= (PINAHGNC[1]!='SUMO3')
FilterArray *= (PINAHGNC[1]!='SUMO4')

# Output (format "ID1\tID2\n"):
Indices = scipy.where(FilterArray)[0]
fw      = open('PINAOnlyHGNCNoSUMONoUBC.tsv','w')
fw.write('NodeA_HGNC\tNodeB_HGNC\n')
for i in Indices:
    fw.write(PINAHGNC[0,i]+'\t'+
             PINAHGNC[1,i]+'\n')
fw.close()

# Filter out a->b;b->a edges:
# Output (format "ID1\tID2\n"):
PINAOnlyHGNCNoSUMONoUBC = scipy.genfromtxt(fname='PINAOnlyHGNCNoSUMONoUBC.tsv',
                                           dtype=str,
                                           delimiter='\t\t', # workaround to have no delimiter :-)
                                           skip_header=1,
                                           unpack=True)
PINAOnlyHGNCNoSUMONoUBC = scipy.sort(PINAOnlyHGNCNoSUMONoUBC)

SortedEdges = []
for i in xrange(len(PINAOnlyHGNCNoSUMONoUBC)):
    SortedEdges.append('\t'.join(sorted(PINAOnlyHGNCNoSUMONoUBC[i].split('\t'))))
SortedEdges = scipy.array(SortedEdges)
UniqEdges,\
UniqIndices = scipy.unique(SortedEdges,
                           return_index=True)

fw      = open('PINAOnlyHGNCNoSUMONoUBCNonRedundant.tsv','w')
fw.write('NodeA_HGNC\tNodeB_HGNC\n')
for i in scipy.sort(UniqIndices):
    fw.write(SortedEdges[i]+'\n')
fw.close()

# Compare Mohammed's HGNC PINA network to Rene's one:
PINAMohammed = scipy.genfromtxt(fname='PINA_noSUMO_noUBC_edges.txt',
                                dtype=str,
                                delimiter='\t\t\t',  # workaround to have no delimiter :-)
                                skip_header=0,
                                unpack=True)

SortedPINAMohammed = []
for i in xrange(len(PINAMohammed)):
    SortedPINAMohammed.append('\t'.join(sorted(PINAMohammed[i].split('\t'))))
SortedPINAMohammed = scipy.sort(scipy.array(SortedPINAMohammed))

PINARene = scipy.genfromtxt(fname='PINAOnlyHGNCNoSUMONoUBCNonRedundant.tsv',
                            dtype=str,
                            delimiter='\t\t\t',  # workaround to have no delimiter :-)
                            skip_header=1,
                            unpack=True)

SortedPINARene = []
for i in xrange(len(PINARene)):
    SortedPINARene.append('\t'.join(sorted(PINARene[i].split('\t'))))
SortedPINARene = scipy.sort(scipy.array(SortedPINARene))

PINAReneGeneSymbolsNotInPINAMohammed = scipy.setdiff1d(ar1=SortedPINARene,
                                                       ar2=SortedPINAMohammed,
                                                       assume_unique=False)
fw = open('PINAReneEntriesNotInPINAMohammed.tsv','w')
for i in xrange(len(PINAReneGeneSymbolsNotInPINAMohammed)):
    fw.write(PINAReneGeneSymbolsNotInPINAMohammed[i]+'\n')
fw.close()

PINAMohammedGeneSymbolsNotInPINARene = scipy.setdiff1d(ar1=SortedPINAMohammed,
                                                       ar2=SortedPINARene,
                                                       assume_unique=False)
fw = open('PINAMohammedEntriesNotInPINARene.tsv','w')
for i in xrange(len(PINAMohammedGeneSymbolsNotInPINARene)):
    fw.write(PINAMohammedGeneSymbolsNotInPINARene[i]+'\n')
fw.close()

SortedPINAOnlyHGNCNoSUMONoUUBC = []
for i in xrange(len(PINAOnlyHGNCNoSUMONoUBC)):
    SortedPINAOnlyHGNCNoSUMONoUUBC.append('\t'.join(sorted(PINAOnlyHGNCNoSUMONoUBC[i].split('\t'))))
SortedPINAOnlyHGNCNoSUMONoUUBC = scipy.sort(scipy.unique(scipy.array(SortedPINAOnlyHGNCNoSUMONoUUBC)))

PINAMohammedGeneSymbolsNotInPINAOnlyHGNCNoSUMONoUBC = scipy.setdiff1d(ar1=SortedPINAMohammed,
                                                                      ar2=SortedPINAOnlyHGNCNoSUMONoUUBC,
                                                                      assume_unique=False)
fw = open('PINAMohammedEntriesNotInPINAOnlyHGNCNoSUMONoUBC.tsv','w')
for i in xrange(len(PINAMohammedGeneSymbolsNotInPINAOnlyHGNCNoSUMONoUBC)):
    fw.write(PINAMohammedGeneSymbolsNotInPINAOnlyHGNCNoSUMONoUBC[i]+'\n')
fw.close()

SortedPINAOnlyHGNC = []
for i in xrange(len(PINAHGNC[0])):
    if(PINAHGNC[0,i]=='UBC'):
        print PINAHGNC[0,i],PINAHGNC[1,i]
    if(PINAHGNC[1,i]=='UBC'):
        print PINAHGNC[0,i],PINAHGNC[1,i]
    SortedPINAOnlyHGNC.append('\t'.join(sorted([PINAHGNC[0,i],PINAHGNC[1,i]])))
SortedPINAOnlyHGNC = scipy.sort(scipy.unique(scipy.array(SortedPINAOnlyHGNC)))

PINAMohammedGeneSymbolsNotInPINAOnlyHGNC = scipy.setdiff1d(ar1=SortedPINAMohammed,
                                                           ar2=SortedPINAOnlyHGNC,
                                                           assume_unique=False)
fw = open('PINAMohammedEntriesNotInPINAOnlyHGNC.tsv','w')
for i in xrange(len(PINAMohammedGeneSymbolsNotInPINAOnlyHGNC)):
    fw.write(PINAMohammedGeneSymbolsNotInPINAOnlyHGNC[i]+'\n')
fw.close()
