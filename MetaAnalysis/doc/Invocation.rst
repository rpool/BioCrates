Invocation of ``QCGWAData.py``
==============================

The module ``QCGWAData.py`` is invoked as follows:

**Usage**::

    usage: QCGWAData.py [-h] [-P NAME]
    
    This python module can be used for the initial QC of GWA output files.
    
    optional arguments:
      -h, --help            show this help message and exit
      -P NAME, --protocolfile NAME
                            NAME: Name of file contains the necessary input and
                            protocol values (input, format XML)

In practice you will run the following command:

**Command**::
   
    ./QCGWAData.py -P <ProtocolFileName>.xml

in your working directory, after linking to the ``QCGWAData.py`` executable.

**Link**::

    ln -s <PathToQCGWADataSource>/QCGWAData.py ./

The ``.xml`` input file
=======================
The module ``QCGWAData.py`` needs an input file in ``.xml`` format that contains all necessary instructions for a successful run: the *protocol*. Below follows a description of the ingredients of this input file. Knowledge on this input format is sufficient for a QC run. For special functionality, you should either modify the ``python`` source code or consult the developer (i.e. René Pool).

An ``.xml`` file consists of tags, specified as **<Tag>** *Value or Child* **</Tag>**. A tag contains values and/or children tags. Each tag has a parent, except for the *root* tag. Module ``QCGWAData.py`` searches for specific tags in the ``.xml`` file. These tags contain the instructions ``QCGWAData.py`` needs for a successful run.

The ``.xml`` needs the following tag blocks (comments are between ``<!-- -->``)

**TagBlocks**::

    <!-- This contains the QC protocol for the TwinsUK GWA output -->
    <Protocol>
        <!-- The directory containing the GWA output files -->
        <GWAOutputPath></GWAOutputPath>

        <!-- Your working directory -->
        <WorkingPath></WorkingPath>

        <!-- Format specifiers of GWA output and GWA extra files -->
        <Format></Format>

        <!-- Hapmap mapping details -->
        <HapMap></HapMap>

        <!-- Files containing extra info associated with the GWA output files -->
        <ExtraInfoFiles>
            <ExtraInfoFile1></ExtraInfoFile1>
            <ExtraInfoFile2></ExtraInfoFile2>
        </ExtraInfoFiles>
        <ExtraInfoColumns>
            <!-- Column specifiers of the extra info file -->
            <!-- ColumnName:        Expected name of the column -->
            <!-- boSetColumnName:   Change the column name? -->
            <!-- SetColumnName:     Change the column name -->
        </ExtraInfoColumns>

        <!-- Details on the ExcludedMetabolites file -->
        <ExcludedMtbFile></ExcludedMtbFile>

        <!-- Details on the PositiveControl file -->
        <PositiveControlFile></PositiveControlFile>
        
        <!-- File containing the list of metabolite names you want to QC -->
        <MtbNameFile></MtbNameFile>
        <MtbGWAColumns>
            <!-- GWA file column specifiers -->
            <!-- ColumnName:            Expected name of the column -->
            <!-- boSetColumnName:       Change the column name? -->
            <!-- SetColumnName:         Change the column name -->
            <!-- MandatoryFieldEntries: Conditions the fields should apply to -->
            <!-- NonNAFieldType:        Data type of the non-NA fields -->
            <!-- Filters:               What filters can be applied? -->
            <!-- NumericFieldPrecision: Precision specifier required for some columns -->
        </MtbGWAColumns>

        <QCChecks>
            <!-- IN: value present in list -->
        </QCChecks>

        <QCFilters>
            <!-- Filter specifiers -->
            <!-- NE: != -->
            <!-- EQ: == -->
            <!-- LT: < -->
            <!-- LE: <= -->
            <!-- GT: > -->
            <!-- GE: >= -->
        </QCFilters>
        <EMAC></EMAC>
    </Protocol>

As we can see the *root* tag should be called ``Protocol``. Its children should be called ``GWAOutputPath``, ``WorkingPath``, ``Format``, ``HapMap``, ``ExtraInfoFiles``, ``ExtraInfoColumns``, ``ExcludedMtbFile``, ``PositiveControlFile``, ``MtbNameFile``, ``MtbGWAColumns``, ``QCChecks``, ``QCFilters`` and ``EMAC``. A brief description for each of the children tags will be given below.

The ``GWAOutputPath`` block
---------------------------
In this block we specify the location of the GWAS data files that need to be QCed.

**GWAOutputPath**::

    <GWAOutputPath>
        /mnt/ntr.biocrates@fluke.psy.vu.nl/PreMAQC/ReneAndHarmen_25Jan2013/TestTwinsUK/Data
    </GWAOutputPath>

The ``WorkingPath`` block
-------------------------
In this block we specify the working directory.

**WorkingPath**::

    <WorkingPath>
        /mnt/ntr.biocrates@fluke.psy.vu.nl/PreMAQC/ReneAndHarmen_25Jan2013/TestTwinsUK/Jobs/Job_0
    </WorkingPath>

The ``Format`` block
--------------------
In this block we specify the format properties of the 'extra info' files and of the GWAS files. Child tag ``ColumnNames`` contains a comma-separated list of what column names to expect, ``Delimiter`` of what delimiter to expect. The latter can be overriden by setting ``boSetDelimeter`` to ``True`` and specifying the new delimeter in ``SetDelimeter``. ``boRemoveDuplicateLines`` specifies whether or not duplicate lines are removed. Child tag ``ExtraInfoDelimiter`` sets the delimiter of the 'extra info' files.

**Format**::

   <Format>
        <ColumnNames>SNPID,chr,position,coded_all,noncoded_all,strand_genome,beta,SE,pval,AF_coded_all,HWE_pval,n_total,imputed,used_for_imp,oevar_imp</ColumnNames>
        <Delimiter>WhiteSpace</Delimiter>
        <boSetDelimeter>False</boSetDelimeter>
        <SetDelimeter>WhiteSpace</SetDelimeter>
        <boRemoveDuplicateLines>True</boRemoveDuplicateLines>
        <ExtraInfoDelimiter>,</ExtraInfoDelimiter>
    </Format> 

The ``HapMap`` block
--------------------
If one wants to map the SNPs listed in a GWAS file to a specific Hapmap build/release, this block can be used to specify the details. ``SummaryPath`` sets the path where your target Hapmap file resides. This file should list the following information: SNPID, chr, position, AlleleA, AlleleB, strand_genome and MAF. ``SummaryFile`` sets the Hapmap summary file name. ``SourceBuild`` specifies the build number of your GWAS file, ``SourceRelease`` its release number and ``SourceRefPanel`` the source reference panel. ``DestBuild`` specifies the target build number of your GWAS file, ``DestRelease`` its release number and ``DestRefPanel`` the target reference panel. ``Delimiter`` specifies the column delimiter in the Hapmap summary file.

**HapMap**::

   <HapMap>
        <SummaryPath>/mnt/ntr.biocrates@fluke.psy.vu.nl/PreMAQC/ReneAndHarmen_25Jan2013/TestTwinsUK/Data</SummaryPath>
        <SummaryFile>HapMapB36R22.txt.gz</SummaryFile>
        <SourceBuild>36</SourceBuild>
        <SourceRelease>22</SourceRelease>
        <SourceRefPanel>CEU</SourceRefPanel>
        <DestBuild>36</DestBuild>
        <DestRelease>22</DestRelease>
        <DestRefPanel>CEU</DestRefPanel>
        <Delimiter>,</Delimiter>
    </HapMap> 

The ``ExtraInfoFiles`` block
----------------------------


    <ExtraInfoFiles>
        <ExtraInfoFile1>
            <boUse>True</boUse>
            <Path>/mnt/ntr.biocrates@fluke.psy.vu.nl/PreMAQC/ReneAndHarmen_25Jan2013/TestTwinsUK/Data</Path>
            <Name>QUAL1-biocrates.csv.gz</Name>
        </ExtraInfoFile1>
        <ExtraInfoFile2>
            <boUse>True</boUse>
            <Path>/mnt/ntr.biocrates@fluke.psy.vu.nl/PreMAQC/ReneAndHarmen_25Jan2013/TestTwinsUK/Data</Path>
            <Name>QUAL2-biocrates.csv.gz</Name>
        </ExtraInfoFile2>
        <ExtraInfoFile3>
            <boUse>True</boUse>
            <Path>/mnt/ntr.biocrates@fluke.psy.vu.nl/PreMAQC/ReneAndHarmen_25Jan2013/TestTwinsUK/Data</Path>
            <Name>QUAL3-biocrates.csv.gz</Name>
        </ExtraInfoFile3>
        <ExtraInfoFile4>
            <boUse>True</boUse>
            <Path>/mnt/ntr.biocrates@fluke.psy.vu.nl/PreMAQC/ReneAndHarmen_25Jan2013/TestTwinsUK/Data</Path>
            <Name>QUAL4-biocrates.csv.gz</Name>
        </ExtraInfoFile4>
    </ExtraInfoFiles>


.. toctree::
   :maxdepth: 2

