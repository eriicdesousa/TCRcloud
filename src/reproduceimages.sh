pip3 install TCRcloud
TCRcloud testdata
TCRcloud download alpharepertoire.airr.json
TCRcloud download betarepertoire.airr.json
TCRcloud cloud alpharepertoire.airr.rearrangements.tsv -l false
TCRcloud cloud betarepertoire.airr.rearrangements.tsv -l false
TCRcloud radar alpharepertoire.airr.rearrangements.tsv -l false -ut 2000 -cht 2000 -ct 0.1
TCRcloud radar betarepertoire.airr.rearrangements.tsv -l false -ut 2000 -cht 2000 -ct 0.1
TCRcloud surface alpharepertoire.airr.rearrangements.tsv -zt 3.5 -yt 25 -yb 5 
TCRcloud surface betarepertoire.airr.rearrangements.tsv -zt 3.5 -yt 28 -yb 8 