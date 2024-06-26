![TCRcloud](https://github.com/eriicdesousa/TCRcloud/raw/main/images/TCRcloud.png)

![GitHub last commit](https://img.shields.io/github/last-commit/eriicdesousa/TCRcloud)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/eriicdesousa/TCRcloud)
![PyPI](https://img.shields.io/pypi/v/TCRcloud)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/TCRcloud)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/TCRcloud)
![License](https://img.shields.io/github/license/eriicdesousa/TCRcloud)

## TCRcloud is an Adaptive Immune Receptor Repertoire (AIRR) visualization and comparison tool

**Instalation**

TCRcloud is written in python and can be installed from PyPI using pip:

    pip3 install TCRcloud

It is compatible with Linux and macOS operating systems and on Windows through the Windows Subsystem for Linux.

TCRcloud uses the AIRR Data Commons API and needs AIRR compliant data as input.

TCRcloud was initially developed for TCR repertoires but it is also compatible with BCR repertoires.

**Usage**

**To create a word cloud**

    TCRcloud cloud repertoire.airr.rearrangements.tsv

By default TCRcloud colours the CDR3 based on the V gene. Only the colours for human TCR and BCR variable genes are coded into TCRcloud but you can provide a json file that atributes colours in Hex format to specific sequences:

    {
    "#FF0000":["CAVSLPTDSWGKLQF","CASSLVVADPYQETQYF"],
    "#0000FF":["CAYRSKGSQGNLIF","CASSLGGQSGNEQFF"]
    }

The sequences not in the json file will be coloured grey.

**To use your custom colours for the word cloud**

    TCRcloud cloud repertoire.airr.rearrangements.tsv -c colours.json

**To create a word cloud without a legend**

    TCRcloud cloud repertoire.airr.rearrangements.tsv -l False

**To create a radar plot comparing diversity indices**

    TCRcloud radar repertoire.airr.rearrangements.tsv

By default TCRcloud uses repertoire_id but you can create a legend with the text you want by providing a json file:

    {
    "PRJNA509910-su008_pre-TRA":"Subject 8 pre-treatment",
    "PRJNA509910-su008_post-TRA":"Subject 8 post-treatment",
    "PRJNA509910-su008_pre-TRB":"Subject 8 pre-treatment",
    "PRJNA509910-su008_post-TRB":"Subject 8 post-treatment"
    }

**To create a radar plot with your desired legend**

	TCRcloud radar repertoire.airr.rearrangements.tsv -c legend.json

**To create a radar plot without a legend**
    
    TCRcloud radar repertoire.airr.rearrangements.tsv -l False

**To export the calculated indices from the radar to a text file**

    TCRcloud radar repertoire.airr.rearrangements.tsv -e True

**The indices utilized in the radar**

**Distinct CDR3:** Count of the unique CDR3 sequences in the sample

**Convergence:** Frequency of CDR3 amino acid sequences that are coded by more than one nucleotide sequence 

**D50 Index:** Developed by Dr Jian Han, this index represents the percent of dominant and unique T or B cell clones that account for the cumulative 50% of the total CDR3 counted in the sample

**Gini Coefficient:** Originally developed by Dr Corrado Gini to represent the wealth inequality within a social group

**Shannon Index:** Originally developed by Dr Claude Shannon to quantify the entropy in strings of text. Calculated here using log<sub>2<sub>

**Gini-Simpson Index:** This is a transformation of the index originally developed by Dr Edward H. Simpson to measure diversity in a ecosystem

**Chao1 index:** Originally developed by Dr Anne Chao to estimate richness in a ecological community

**To create an amino V-genes plot**

    TCRcloud vgenes repertoire.airr.rearrangements.tsv

**To export the processed data from the V-genes plot to a text file**

    TCRcloud vgenes repertoire.airr.rearrangements.tsv -e True

**To create a 2-D amino acid plot**

    TCRcloud aminoacids repertoire.airr.rearrangements.tsv

**To create a 3-D amino acid plot**

    TCRcloud aminoacids repertoire.airr.rearrangements.tsv -t True

**To export the processed data from the amino acid plot to a text file**

    TCRcloud aminoacids repertoire.airr.rearrangements.tsv -e True

Using TCRcloud you can download rearragements files from the AIRR compliant databases based on AIRR repertoire metadata files

**To download AIRR rearrangements files**

    TCRcloud download repertoire.airr.json

TCRcloud provides some test data to experiment the tool. The data is from subject number 8 of Yost et al. (2019) (DOI:[ 10.1038/s41591-019-0522-3](https://doi.org/10.1038/s41591-019-0522-3))

**To download the test data repertoire file**

    TCRcloud testdata

After having the alpharepertoire.airr.json  and betarepertoire.airr.json file you can use the download function included in TCRcloud to get the matching rearragements file.

## Examples:

**TRA CDR3 word cloud**
![alpha](https://github.com/eriicdesousa/TCRcloud/raw/main/images/alpha.png)

**TRB CDR3 word cloud**
![beta](https://github.com/eriicdesousa/TCRcloud/raw/main/images/beta.png)

**Diversity comparison**
![radar](https://github.com/eriicdesousa/TCRcloud/raw/main/images/radar.png)

**V-genes plot**
![vgene](https://github.com/eriicdesousa/TCRcloud/raw/main/images/vgene.png)

**2-D Amino acid plot**
![2amino](https://github.com/eriicdesousa/TCRcloud/raw/main/images/2amino.png)

**3-D Amino acid plot**
![3amino](https://github.com/eriicdesousa/TCRcloud/raw/main/images/3amino.png)