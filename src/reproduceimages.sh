pip3 install TCRcloud
TCRcloud testdata
TCRcloud download alpharepertoire.airr.json
TCRcloud download betarepertoire.airr.json
TCRcloud cloud alpharepertoire.airr.rearrangements.tsv -l false
TCRcloud cloud betarepertoire.airr.rearrangements.tsv -l false
TCRcloud radar alpharepertoire.airr.rearrangements.tsv -l false
TCRcloud radar betarepertoire.airr.rearrangements.tsv -l false
TCRcloud surface -zt 3.5 -yt 25 -yb 5 alpharepertoire.airr.rearrangements.tsv
TCRcloud surface -zt 3.5 -yt 28 -yb 8 betarepertoire.airr.rearrangements.tsv