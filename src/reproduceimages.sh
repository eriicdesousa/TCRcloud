pip3 install TCRcloud
TCRcloud testdata
TCRcloud download alpharepertoire.airr.json
TCRcloud download betarepertoire.airr.json
TCRcloud cloud alpharepertoire.airr.rearrangements.tsv
TCRcloud cloud betarepertoire.airr.rearrangements.tsv
TCRcloud radar alpharepertoire.airr.rearrangements.tsv -l false -uh 2000 -chh 2000 -ch 0.1
TCRcloud radar betarepertoire.airr.rearrangements.tsv -l false -uh 2000 -chh 2000 -ch 0.1
TCRcloud vgenes alpharepertoire.airr.rearrangements.tsv -zla -0.01 -zha 4 -yha 23 -yla 5 
TCRcloud vgenes betarepertoire.airr.rearrangements.tsv -zlb -0.01 -zhb 4 -yhb 28
TCRcloud vgenes alpharepertoire.airr.rearrangements.tsv -c True
TCRcloud vgenes betarepertoire.airr.rearrangements.tsv -c True
TCRcloud aminoacids alpharepertoire.airr.rearrangements.tsv -l 21
TCRcloud aminoacids betarepertoire.airr.rearrangements.tsv -l 26
TCRcloud aminoacids alpharepertoire.airr.rearrangements.tsv -l 21 -t True
TCRcloud aminoacids betarepertoire.airr.rearrangements.tsv -l 26 -t True
TCRcloud aminoacids alpharepertoire.airr.rearrangements.tsv -l 21 -t True -c True
TCRcloud aminoacids betarepertoire.airr.rearrangements.tsv -l 26 -t True -c True