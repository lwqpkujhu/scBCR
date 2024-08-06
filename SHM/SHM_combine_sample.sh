# # parse column header
# awk 'NR == 1 {print $0}' ./21_Tube_3_results/heavy_parse-select.tsv > combined_heavy_parse-select.tsv
# 
# # parse all samples to one file
# for i in `ls -d */`; do 
# awk 'NR >= 2 {print $0}' ./${i}heavy_parse-select.tsv >> combined_heavy_parse-select.tsv;
# done


# DefineClones
DefineClones.py -d combined_heavy_parse-select.tsv --act set --model ham \
    --norm len --dist 0.16

# germlines by clone
CreateGermlines.py -d combined_heavy_parse-select_clone-pass.tsv -g dmask --cloned \
    -r /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHD.fasta /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta


# germlines NOT by clone
CreateGermlines.py -d combined_heavy_parse-select_clone-pass.tsv -g dmask \
    -r /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHD.fasta /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta
