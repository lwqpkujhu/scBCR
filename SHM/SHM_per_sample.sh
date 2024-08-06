## process each sample
AssignGenes.py igblast -s ../$1/$1/vdj_b/filtered_contig.fasta -b /home/weiqiang/software/ncbi-igblast-1.14.0/share/igblast \
   --organism mouse --loci ig --format blast --outdir ./$1

MakeDb.py igblast -i ./$1/filtered_contig_igblast.fmt7 -s ../$1/$1/vdj_b/filtered_contig.fasta \
   -r /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt/mouse/vdj/imgt_mouse_*.fasta \
   --10x ../$1/$1/vdj_b/filtered_contig_annotations.csv --extended

cd $1/

ParseDb.py select -d filtered_contig_igblast_db-pass.tsv -f locus -u "IGH" \
        --logic all --regex --outname heavy
ParseDb.py select -d filtered_contig_igblast_db-pass.tsv -f locus -u "IG[LK]" \
        --logic all --regex --outname light

DefineClones.py -d heavy_parse-select.tsv --act set --model ham \
    --norm len --dist 0.16

CreateGermlines.py -d heavy_parse-select_clone-pass.tsv -g dmask --cloned \
    -r /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHD.fasta /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta /home/weiqiang/software/ncbi-igblast-1.14.0/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta



