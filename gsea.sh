#!/bin/bash
me=whoami

### ============= no need to modify
folder=/public/$me/projects/CMLblast
gsea_path=/public/$me/utilities/gsea-4.4.0.jar

cd $folder
mdkir $OUTDIR

java -cp $gsea_path \
    -Xmx5g xtools.gsea.GseaPreranked \
    -rpt_label $LABEL \
    -rnk $folder/$RNKFILE \
    -gmx $GMT \
    -out $folder/$OUTDIR \
    -plot_top_x 100 \
    -collapse false \
    -mode Abs_max_of_probes \
    -norm meandiv \
    -scoring_scheme weighted \
    -include_only_symbols true \
    -make_sets true \
    -rnd_seed 149 \
    -zip_report false \
    -gui false \
    -nperm 10000 \
    -set_min 15 \
    -set_max 500