![TCRcloud](https://github.com/oldguyeric/TCRcloud/raw/main/images/TCRcloud.png)

## TCRcloud is a TCR repertoire visualization and comparison tool

**Instalation**

TCRcloud is written in python and can be installed from PyPI using pip

    $ pip3 install TCRcloud

TCRcloud uses the AIRR Data Commons API and needs AIRR compliant data as input

**Usage**

**To create a wordcloud**

    $ TCRcloud cloud -a AIRR_rearrangements_file.tsv

By default TCRcloud colours the CDR3 based on the Vgene but you can provide a json file that atributes colours in Hex format to specific sequences:

    {
    "#FF0000":["CAASITGNQFYF","CAVREDGTSGSARQLTF"],
    "#0000FF":["CAVMDSNYQLIW"]
    }

The sequences not in the file will be coloured grey.

**To use your custom colours for the word cloud**

    $ TCRcloud cloud -a AIRR_rearrangements_file.tsv -c colours.json

**To create a radar plot comparing diversity metrics**

    $ TCRcloud radar -a AIRR_rearrangements_file.tsv

By default TCRcloud uses repertoire_id but you can create a legend with the text you want by providing a json file:

    {
    "2839362682105696746-242ac113-0001-012":"Twin 2A",
    "2939134772391776746-242ac113-0001-012":"Twin 2B"
    }  

**To create a radar plot with your desired legend**

	$ TCRcloud radar -a AIRR_rearrangements_file.tsv -l legend.json

Using TCRcloud you can download rearragements files from the AIRR compliant databases based on AIRR repertoire metadata files

**To download AIRR rearrangements files**

    $ TCRcloud download -a AIRR_repertoire_file.json

TCRcloud provides some test data to experiment the tool. The data is one twin pair from the monozygotic twins study from the Mark Davis lab (DOI:[ 10.1038/ncomms11112](https://doi.org/10.1038/ncomms11112))

**To download the test data repertoire file**

    $ TCRcloud testdata

After having the testdata.airr.json file you can use the download function included in TCRcloud to get the matching rearragements file.

## Examples:

**TRA CDR3 word cloud**
![alpha](https://github.com/oldguyeric/TCRcloud/raw/main/images/alpha.png)

**TRB CDR3 word cloud**
![beta](https://github.com/oldguyeric/TCRcloud/raw/main/images/beta.png)

**TRG CDR3 word cloud**
![gamma](https://github.com/oldguyeric/TCRcloud/raw/main/images/gamma.png) 

**TRD CDR3 word cloud**
![delta](https://github.com/oldguyeric/TCRcloud/raw/main/images/delta.png) 

**Diversity comparison**

Comparing here the αβ repertoire from one twin pair from the monozygotic twins study from the Mark Davis lab (DOI:[ 10.1038/ncomms11112](https://doi.org/10.1038/ncomms11112))

![radar](https://github.com/oldguyeric/TCRcloud/raw/main/images/radar.png)
