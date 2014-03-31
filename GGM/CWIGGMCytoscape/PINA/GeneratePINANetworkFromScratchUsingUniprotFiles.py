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

# Read in the mapping table downloaded from www.uniprot.org by querying the all unique UniprotKB Accessions (file "UniprotAC_or_ID_to_EntrezGeneID.tsv")
fr                      = open('UniprotAC_or_ID_to_EntrezGeneID.tsv','r')
UniprotKBID2EntrezIDHdr = fr.readline()
UniprotKBID2EntrezIDHdr = fr.readline()
UniprotKBID2EntrezIDHdr = fr.readline().strip().split('\t')
fr.close()
UniprotKBID2EntrezID    = scipy.genfromtxt(fname='UniprotAC_or_ID_to_EntrezGeneID.tsv',
                                           dtype=str,
                                           delimiter='\t',
                                           skip_header=3,
                                           unpack=True)

# Read in the mapping table of EntrezID to HGNC symbol downloaded from ensembl BioMart (file "BioMartGeneSymbol2EntrezID.tsv")
fr               = open('BioMartGeneSymbol2EntrezID.tsv','r')
HGNC2EntrezIDHdr = fr.readline().strip().split('\t')
fr.close()
HGNC2EntrezID    = scipy.genfromtxt(fname='BioMartGeneSymbol2EntrezID.tsv',
                                    dtype=str,
                                    delimiter='\t',
                                    skip_header=1,
                                    unpack=True)
# Write translated UniprotKB IDs:
if(False):
    fw = open('UsingUniprotFiles/PINAHGNC.tsv','w')
    fw.write('NodeA_HGNC\tNodeB_HGNC\tNodeAUniprotKB\tNodeBUniprotKB\n')
    for i in xrange(len(PINAUniprot[0])):
        NodeA          = PINAUniprot[0,i]
        NodeB          = PINAUniprot[1,i]
        NodeAEntrezIDs = UniprotKBID2EntrezID[1,scipy.where(UniprotKBID2EntrezID[0]==NodeA)[0]]
        NodeAEntrezIDs = NodeAEntrezIDs[scipy.where(NodeAEntrezIDs!='')[0]]
        NodeBEntrezIDs = UniprotKBID2EntrezID[1,scipy.where(UniprotKBID2EntrezID[0]==NodeB)[0]]
        NodeBEntrezIDs = NodeBEntrezIDs[scipy.where(NodeBEntrezIDs!='')[0]]
        NodeAHGNCSmbls = ['None']
        for j in xrange(len(NodeAEntrezIDs)):
            EntrezID       = NodeAEntrezIDs[j]
            HGNCSmbls      = scipy.unique(HGNC2EntrezID[1,scipy.where(HGNC2EntrezID[0]==EntrezID)[0]])
            HGNCSmbls      = HGNCSmbls[scipy.where(HGNCSmbls!='')[0]]
            NodeAHGNCSmbls = HGNCSmbls.tolist()
        NodeBHGNCSmbls = ['None']
        for j in xrange(len(NodeBEntrezIDs)):
            EntrezID       = NodeBEntrezIDs[j]
            HGNCSmbls      = scipy.unique(HGNC2EntrezID[1,scipy.where(HGNC2EntrezID[0]==EntrezID)[0]])
            HGNCSmbls      = HGNCSmbls[scipy.where(HGNCSmbls!='')[0]]
            NodeBHGNCSmbls = HGNCSmbls.tolist()
        Product = list(itertools.product(NodeAHGNCSmbls,NodeBHGNCSmbls))
        for j in xrange(len(Product)):
            fw.write(Product[j][0]+'\t'+
                     Product[j][1]+'\t'+
                     NodeA+'\t'+
                     NodeB+'\t\n')
    fw.close()

# Find unmatched UniprotKBIDs
PINAHGNC = scipy.genfromtxt(fname='UsingUniprotFiles/PINAHGNC.tsv',
                            dtype=str,
                            delimiter='\t',
                            skip_header=1,
                            unpack=True)
UnmatchedUniprotKBIDs = scipy.array([])
UnmatchedUniprotKBIDs = scipy.append(UnmatchedUniprotKBIDs,PINAHGNC[2,scipy.where(PINAHGNC[0]=='None')[0]])
UnmatchedUniprotKBIDs = scipy.append(UnmatchedUniprotKBIDs,PINAHGNC[3,scipy.where(PINAHGNC[1]=='None')[0]])
UnmatchedUniprotKBIDs = scipy.unique(UnmatchedUniprotKBIDs)
fw = open('UsingUniprotFiles/UnmatchedUniprotKBIDs.txt','w')
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
fw = open('UsingUniprotFiles/ENGAGEGeneSymbolsNotInPINA.txt','w')
for i in xrange(len(ENGAGEGeneSymbolsNotInPINA)):
    fw.write(ENGAGEGeneSymbolsNotInPINA[i]+'\n')
fw.close()

PINAGeneSymbolsNotInENGAGE = scipy.setdiff1d(ar1=GeneSymbolsInPINA,
                                             ar2=UniqHGNCSymbolsInENGAGEData,
                                             assume_unique=True)
fw = open('UsingUniprotFiles/PINAGeneSymbolsNotInENGAGE.txt','w')
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
fw      = open('UsingUniprotFiles/PINAOnlyHGNCNoSUMONoUBC.tsv','w')
fw.write('NodeA_HGNC\tNodeB_HGNC\n')
for i in Indices:
    fw.write(PINAHGNC[0,i]+'\t'+
             PINAHGNC[1,i]+'\n')
fw.close()

# Filter out a->b;b->a edges:
# Output (format "ID1\tID2\n"):
PINAOnlyHGNCNoSUMONoUBC = scipy.genfromtxt(fname='UsingUniprotFiles/PINAOnlyHGNCNoSUMONoUBC.tsv',
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

fw      = open('UsingUniprotFiles/PINAOnlyHGNCNoSUMONoUBCNonRedundant.tsv','w')
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

PINARene = scipy.genfromtxt(fname='UsingUniprotFiles/PINAOnlyHGNCNoSUMONoUBCNonRedundant.tsv',
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
fw = open('UsingUniprotFiles/PINAReneEntriesNotInPINAMohammed.tsv','w')
for i in xrange(len(PINAReneGeneSymbolsNotInPINAMohammed)):
    fw.write(PINAReneGeneSymbolsNotInPINAMohammed[i]+'\n')
fw.close()

PINAMohammedGeneSymbolsNotInPINARene = scipy.setdiff1d(ar1=SortedPINAMohammed,
                                                       ar2=SortedPINARene,
                                                       assume_unique=False)
fw = open('UsingUniprotFiles/PINAMohammedEntriesNotInPINARene.tsv','w')
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
fw = open('UsingUniprotFiles/PINAMohammedEntriesNotInPINAOnlyHGNCNoSUMONoUBC.tsv','w')
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
fw = open('UsingUniprotFiles/PINAMohammedEntriesNotInPINAOnlyHGNC.tsv','w')
for i in xrange(len(PINAMohammedGeneSymbolsNotInPINAOnlyHGNC)):
    fw.write(PINAMohammedGeneSymbolsNotInPINAOnlyHGNC[i]+'\n')
fw.close()