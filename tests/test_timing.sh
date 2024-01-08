#!/bin/bash

# this is just a template to test timing of the subcommands

# nontr
python -m cProfile -s 'time' -m snoopsv nontr -v $VCF_IN -o $VCF_OUT -s $SAMPLE -b $BAM -c chr14 -i 0 -n 10 | tee out_nontr.txt

# tr
python -m cProfile -s 'time' -m snoopsv tr -a $ANNOT -v $VCF_IN -o $VCF_OUT -s $SAMPLE -b $BAM --exclude-columns SEQ --annot-columns CHROM,START,End,SRC,REP_LEN,REP_CNT,SEQ -c chr1 -i 0 -n 100 | tee out_tr.txt
