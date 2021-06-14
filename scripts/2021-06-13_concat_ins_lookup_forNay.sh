# prepare list of instruments for lookup

# instrument file has been previously prepared from BMI-HF project (see Scratch folder)
FILE_INS="data/BMI_MVMR_instruments_rsq0.01/instrument_chr_bp.tsv"
FILE_GWAS="data/GWAS/BMI_GIANT+UKB2018.CLEAN.tsv.gz"
OUTFILE="data/BMI-AF-CAD-T2D-SBP_instrument_rsq0.01.tsv"

# after tabix lookup 1704 of 1706 instrument data are available in the GWAS file
tabix -R $FILE_INS $FILE_GWAS | \
  awk 'BEGIN {OFS="\t"; print "CHR", "BP", "rsID", "EA", "OA"}
       {print $1, $2, $3, $4, $5}' > $OUTFILE
