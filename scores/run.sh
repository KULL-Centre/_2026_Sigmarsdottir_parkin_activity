
wget https://sid.erda.dk/share_redirect/ff6CDV5565/facs_data.tgz
tar -zxf facs_data.tgz
Rscript facs.r facs_data/*.csv

Rscript scores.r ../counts/counts_prkn_bc.csv > scores.out

