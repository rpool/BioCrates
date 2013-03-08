import argparse

import Logger

def LogArguments(Log=Logger,
                 ArgParser=argparse.ArgumentParser,
                 Arguments=argparse.Namespace):
    ArgParser.print_help()
    ArgParser.print_help(Log.GetFileHandle())

#   Calculate max. of keylength for formatting
    MaxLen = 0
    for Key in vars(Arguments).iterkeys():
        MaxLen = max(MaxLen,len(Key))
    LogString =  '\n****\n'+\
                 'Used arguments:\n'+\
                 '---------------'
    print LogString
    Log.Write(LogString+'\n')

    FormatString = '{0:<'+str(MaxLen)+'}'
    LogString    = ''
    for Key,Value in vars(Arguments).iteritems():
        LogString += FormatString.format(Key)+': '+str(Value)+'\n'
    LogString += '****\n'
    print LogString
    Log.Write(LogString+'\n')
    return

def ParseArguments(Log=None):
    ArgumentParser = argparse.ArgumentParser(description=\
                                             'This python module can be used for QC of BioCrates metabolomics data.')
    ArgumentParser.add_argument('-e',
                                '--excelfilesample',
                                dest='ExcelSampleDataFileName',
                                help='FILENAME: Name of the excel file containing the sample data (input, \".xls\")',
                                metavar='FILENAME')
    ArgumentParser.add_argument('-E',
                                '--excelfileref',
                                dest='ExcelReferenceDataFileName',
                                help='FILENAME: Name of the excel file containing the reference data (input, \".xls\")',
                                metavar='FILENAME')
    ArgumentParser.add_argument('-s',
                                '--excelsheetnamesample',
                                dest='ExcelSampleSheetName',
                                help='PROPERTY: Name of the excel sheet to be used for the sample data (input)',
                                metavar='PROPERTY')
    ArgumentParser.add_argument('-S',
                                '--excelsheetnameref',
                                dest='ExcelReferenceSheetName',
                                help='PROPERTY: Name of the excel sheet to be used for the reference data (input)',
                                metavar='PROPERTY')
    ArgumentParser.add_argument('-d',
                                '--dlhandling',
                                dest='DLHandling',
                                type=str,
                                help='PROPERTY: Determines the way concentrations below the detection/quantification limits is handled (input).\n'+\
                                    '\"LOD\" sets every concentration below the LOD to \"missing\"\n'+\
                                    '\"LOD-LLOQ\" causes exclusion of a metabolite for which the mean concentration is below LLOQ.',
                                metavar='PROPERTY',
                                default='LOD',
                                choices=['LOD','LOD-LLOQ','OFF'])
    ArgumentParser.add_argument('-C',
                                '--removecustommetabolites',
                                dest='boRemoveCustomMetabolites',
                                help='FLAG: Remove custom metabolites from dataset? (input)',
                                action='store_true',
                                default=False)
    ArgumentParser.add_argument('-p',
                                '--plotbasename',
                                dest='PlotBaseName',
                                help='FILENAME: Base name of the output plots (input,\".pdf\")',
                                default='Plot.pdf',
                                metavar='FILENAME')
    ArgumentParser.add_argument('-o',
                                '--summaryplotname',
                                dest='SummaryPlotBaseName',
                                help='FILENAME: Base name of the file summarizing the plots (input)',
                                default='Plots.pdf',
                                metavar='FILENAME')
    ArgumentParser.add_argument('-a',
                                '--anranges',
                                dest='AnalyticalRangesExcelFileName',
                                help='FILENAME: Name of the excel file containing the analytical ranges of each metabolite (input,\".csv\")',
                                default='SomeFile.csv',
                                metavar='FILENAME')
    ArgumentParser.add_argument('-q',
                                '--QCreference',
                                dest='boQCReference',
                                help='FLAG: Perform QC on reference data (input)',
                                action='store_true',
                                default=False)
    ArgumentParser.add_argument('-Q',
                                '--QCsample',
                                dest='boQCSample',
                                help='FLAG: Perform QC on sample data (input)',
                                action='store_true',
                                default=False)
    ArgumentParser.add_argument('-P',
                                '--plot',
                                dest='boPlotDistributions',
                                help='FLAG: Plot distributions of sampled concentrations (input)',
                                action='store_true',
                                default=False)
    ArgumentParser.add_argument('-i',
                                '--impute',
                                dest='boImpute',
                                help='FLAG: Impute missing data (input)',
                                action='store_true',
                                default=False)
    ArgumentParser.add_argument('-A',
                                '--assumeqcfinished',
                                dest='boQCDone',
                                help='FLAG: Assume that QC and imputation is already performed (input)',
                                action='store_true',
                                default=False)
    ArgumentParser.add_argument('-R',
                                '--correlatecombinations',
                                dest='boCorrelateCombinations',
                                help='FLAG: Determine correlation coefficients of combinations (input)',
                                action='store_true',
                                default=False)
    ArgumentParser.add_argument('-O',
                                '--readrobjectfromfile',
                                dest='boReadRObject',
                                help='FLAG: Read the imputed R-object from file (input)',
                                action='store_true',
                                default=False)
    ArgumentParser.add_argument('-N',
                                '--normalitytest',
                                dest='boNormalityTest',
                                help='FLAG: Perform normality tests on the concentrations, the ln-transformed concentrations and all possible ratios in ln-space (input)',
                                action='store_true',
                                default=False)
    ArgumentParser.add_argument('-z',
                                '--imputezeros',
                                dest='boImputeZeros',
                                help='FLAG: Set zero data entries to \"NA\" in order to impute these values that will become \"missing\" after ln-transformation (input)',
                                action='store_true',
                                default=False)


    Arguments = ArgumentParser.parse_args()

    return ArgumentParser,\
           Arguments