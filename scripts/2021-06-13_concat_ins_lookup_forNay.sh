# prepare list of instruments for lookup

# instrument file has been previously prepared from BMI-HF project (see Scratch folder)
FILE_INS="data/BMI_MVMR_instruments_rsq0.01/Clumping_MVMR/BMI-SBP_instrument"
FILE_GWAS="data/GWAS/BMI_GIANT+UKB2018.CLEAN.tsv.gz"
OUTFILE="data/BMI-SBP_instrument_rsq0.01.tsv"

# after tabix lookup 1704 of 1706 instrument data are available in the GWAS file
tabix -R <(awk -v OFS="\t" 'NR > 1 {print $1, $3}' $FILE_INS | sort -nk2,2 -nk1,1) $FILE_GWAS | \
  awk 'BEGIN {OFS="\t"; print "CHR", "BP", "rsID", "EA", "OA"}
       {print $1, $2, $3, $4, $5}' > $OUTFILE
